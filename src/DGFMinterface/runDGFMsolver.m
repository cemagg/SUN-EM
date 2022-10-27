function [dgfm] = runDGFMsolver(Const, Solver_setup, zMatrices, yVectors, xVectors, mbfs, ngf)
    %runDGFMsolver
    %   Usage:
    %       [dgfm] = runDGFMsolver(Const, Solver_setup, zMatrices, yVectors, xVectors, mbfs, ngf)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details        
    %       zMatrices
    %           The Z-matrices data (can be calculated with e.g. FEKO, or SUN-EM)
    %       yVectors
    %           The Yrhs-vector data
    %       xVectors
    %           The Xsol-vector data (e.g. MoM solution of FEKO or SUN-EM)
    %       mbfs
    %           The MBFs (prim & sec.) generated in runMBFgenerator
    %       ngf
    %           The Numerical Green's Function (NGF) for the static
    %           domain (if used)
    %
    %   Output Arguments:
    %       dgfm
    %           Structs containing DGFM solution and timing data
    %
    %   Description:
    %       Runs the Domain Green's Function Method (DGFM) solution based
    %       on the Z and Y data that was read / parsed from the FEKO *.out,
    %       *.mat, *.str and *.rhs files and also the array data that was setup.
    %       See also FEKDDM-6.1: We need to add a correction term for the
    %       NGF-enhanced DGFM. Also the array forms the dynamic domain and we
    %       need to store a sparse representation of the active impedance matrices.
    %
    %   References:
    %   [1] D. J. Ludick, R. Maaskant, D.B. Davidson, U. Jakobus, R. Mittra and D. de Villiers
    %       "Efficient Analysis of Large Aperiodic Antenna Arrays Using the Domain Green's Function
    %        Method," IEEE Trans. Antennas and Propagation, vol. 62, no. 4, pp. 1-11, 2014.

    narginchk(7,7);
    
    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running DGFM solver'));
    message_fc(Const,sprintf('  (Alpha-weighting scheme: %d and Solution method: %d)', ...
        Const.DGFMweightVectorCalcScheme,Const.useDGFMmethod));
    if (Const.useACA && Const.ACAalg)
        message_fc(Const,sprintf('  (Using the ACA - algorithm: %d)',Const.ACAalg));
    end%if
    if (Const.useDGFMinterpolation)
        message_fc(Const,sprintf('  (Using Interpolation - algorithm: %d with sampling factor %.2f)',...
            Const.useDGFMinterpolation, Const.DGFMinterpolationSamplingFactor));
    end%if
    message_fc(Const,sprintf('  Number of DGFM basis functions: %d.',Solver_setup.num_mom_basis_functions));

    % Initialisations
    dgfm  = [];
    dgfm.name = 'dgfm';
    
    Nmom = Solver_setup.num_mom_basis_functions;  % Total number of basis functions for whole problem
    Nngf = Solver_setup.num_ngf_basis_functions;  % Number of basis functions for NGF domain
    Narr = Nmom;                                  % Total number of basis functions for the array (only as Nngf = 0 for now)
    Ndgfm = Solver_setup.mom_basis_functions_per_array_element;   % Number of basis functions per array element
    numArrayEls = Solver_setup.num_finite_array_elements;  % The number of array elements

    numSols = 1;                                  % The number of reference solutions (TO-DO: only have 1 for now)
    dgfm.numSols = numSols;                       % Calculate a solution for each configuration
    dgfm.Isol = complex(zeros(Narr,numSols));

    % We need to run the DGFM (or i-DGFM) for each frequency point.
    numFreq = Solver_setup.frequencies.freq_num;  % The number of frequency points to process
    numRHSperFreq = 1;                            % The number of solutions per frequency point
                                           
    % Some info about the solution configurations
    message_fc(Const,sprintf('  numSols : %d', numSols));
    message_fc(Const,sprintf('  numFreq : %d', numFreq));
    message_fc(Const,sprintf('  numRHSperFreq : %d', numRHSperFreq));

    % Correct the timing here (now we work with timing arrays) before starting the frequency loop
    dgfm.ngfdgfm_init_time = zeros(1,numFreq);
    dgfm.solTime = zeros(1,numFreq);
    dgfm.totsolTime = 0.0;
    dgfm.totngfdgfm_init_time = 0.0;

    % Calculate the DGFM weighting vectors (note, this is done OUTSIDE the frequency + sol. configuration \
    % loop)
    [dgfm.weightVectors, jackDGFM] = calcDGFMweightVectors(Const, Solver_setup, zMatrices, yVectors, xVectors, mbfs);

    % If we are going to do interpolation afterwards, then we are not going
    % to calculate a full-wave DGFM solution (+ MBF extraction) for each of
    % the elements. Instead, we are going to calculate only a few based on
    % a sampling grid, which we calculate now here:
    
    % TO-DO: Might have to time this initial setup for the interpolation
    if (Const.useDGFMinterpolation)
        % First extract the known array positions
        dgfm.array_XY = zeros(numArrayEls,2);    
        [dgfm.array_XY(:,1),dgfm.array_XY(:,2)] = extract_array_element_positions(Const);
    
        % Now extract some uniformly distributed samples on this grid
        [dgfm.interpolation_sampling_array_indices] = extractDGFMcircularInterpolationGrid(Const, numArrayEls, dgfm.array_XY);
    
        % After we have set up the lattice, we also need to define the reference MBFs:
        % Allocate some space for the MBF coefficients associated with each array element. We need to keep the number of 
        % MBFs per domain consistent. As an initial test, we keep the MBFs for the first element in the array 
        % (which for a circular layout, corresponds to an element near the centre)
        ref_domain = 1;
        dgfm.numRefMBFs = mbfs.numRedMBFs(ref_domain,1);
        ref_domain_basis_functions = Solver_setup.rwg_basis_functions_domains{ref_domain};
        % Note: Below, we assume that we are working with solution index 1
        dgfm.refMBFs = mbfs.RedIsol(ref_domain_basis_functions,1:dgfm.numRefMBFs,ref_domain,1);

        % Initialise complex array for storing the MBF weights
        dgfm.MBF_weights = complex(zeros(numArrayEls,dgfm.numRefMBFs));
    
    end%if
    
    % --------------------- START HERE FREQUENCY LOOP
    for freq=1:numFreq

        % Extract the solutions from which to start and end our DGFM process
        solStart = 1;
        solEnd   = numRHSperFreq;

        % CMA related variables
        % See WORKBOOK 2018-1, p.4 (2018-01-03): Disable now the hybrid CMA
        % DGFM solution here. We now export the Zact impedance matrices and
        % then run the standard CMA part in runEMsolvers.m
%         if (Const.runCMAfromDGFMsolver)
%             % The following flag is used to read the correct *.mat and *.str file
%             Const.cmaActive = true;
%             cma = [];
%             % Read the FEKO reference file in which the CMA results will be returned
%             % for each element individually.
%             % TO-DO: Danie, check this part again!
%             % 2017-12-07: This part has to be read outside of the frequency loop
%             [xVectorsCMA_array_element] = readFEKOXvectorFromFile(Const);
% 
%             % Keep track of the array element IDs, as we append these IDs to the CMA output
%             % files.
%             Const.arrayID = 0;
%         end%if

        % Start timing for NGF-en. DGFM time
        tic
        % See issue FEKDDM-6.1: For the NGF-en. DGFM solver we need to calculate the global coupling
        % modified Zsd' matrix, i.e. Zsd' = inv(Zss)Zsd, where Zsd represents all the coupling matrices
        % between the static and dynamic domains and also the correction
        % matrix, that is Zds_glob*Zsd_glob, see Worbook, pp. 109. Note, the
        % NGF is the same for each solution configuration (except, of course,
        % for frequency solutions)
        % TO-DO: Danie, check here the effect on frequency domain solutions.
        if (Const.runNGFenDGFMsolver)
            Zsd = zMatrices.values(1:Nngf,Nngf+1:Nmom, freq);
            % Back-wards substitution using the NGF already factor L and U components
            b = ngf.L(:,:,freq)\Zsd;
            Zsd_glob = (ngf.U(:,:,freq)\b);
            Zds_glob = zMatrices.values(Nngf+1:Nmom,1:Nngf, freq);
            %Zcorr = Zms*Zsd_glob(:,(m-1)*Ndgfm+1:m*Ndgfm);
            %Zcorr = Zds_glob*Zsd_glob; % <-- size(Narr, Narr)
        end
        dgfm.ngfdgfm_init_time(freq) = toc;

        % Start time again
        tic

        % Allocate space for the local active impedance matrix
        Zact = complex(zeros(Ndgfm,Ndgfm));
        Umn = []; % For the ACA
        Vmn = [];

        % See issue FEKDDM-10: We added now support for multiple solution
        % configurations. Each of the solution configurations has given rise to
        % a suitable weightVectors (based on the excitation-law, exact current,
        % etc. for that specific solution). We need to repeat the DGFM active
        % impedance matrix calculation of each array element for each solution
        % configuration.
        for solNum = solStart:solEnd
            % For a frequency loop, we need to calculate the correct offset for
            % accessing the weighting coefficient.
            index = solNum + (freq-1)*numRHSperFreq;
            for m=1:numArrayEls % M_LOOP
                
                % If we are going to use interpolation, then we only
                % calculate the DGFM entry if it is included in the
                % sampling point
                if (Const.useDGFMinterpolation)
                    if find(dgfm.interpolation_sampling_array_indices == m)
                        % Calculate the DGFM solution for this interpolation sample point.
                        message_fc(Const,sprintf('    Array index %d in sampling grid for interpolation',m));
                    else
                        % Skip this element (solution will be calculated using interpolation)
                        message_fc(Const,sprintf('    Skipping array index %d (not in sampling grid for interpolation)',m));
                        continue;
                    end
                end%if

                % See issue FEKDDM-6.2: Improve the array numbering by using a
                % bottom and top basis function index (particularly for use with
                % the NGF-en. DGFM) - also for the IFB-DGFM (where the finite array is offset in the
                % matrix)
                %domain_bot_m = Const.arrayElBasisFunctRange(m,1);
                %domain_top_m = Const.arrayElBasisFunctRange(m,2);

                domain_m_basis_functions = Solver_setup.rwg_basis_functions_domains{m};

                % A: Calculate the current on the active element if it is active
                % B: Extract the current from the exact solution / MBFs / set
                %    to zero if the element is not active
                if (Const.is_array_element_active(m,solNum))
                    % Solve the array element current by formulating an active
                    % impedance matrix for domain/array element m
                    if (Const.runNGFenDGFMsolver)
                        % For the NGF-en. DGFM, we have here the coupling between the
                        % static domain and the array element being analysed
                        % TO-DO: Check here the effect of a frequency loop.
                        %Zms = zMatrices.values(domain_bot_m:domain_top_m,1:Nngf, freq);
                        Zms = zMatrices.values(domain_m_basis_functions,1:Nngf, freq);
                    end
                    % TO-DO: - Assumed here is 1-to-1 mapping, i.e. Const.arrayMappingVector
                    %          not yet used (would be needed as we consider mixed basis
                    %          functions

                    % The ACA compression of the self-term here cause and error in buildMoMblockACA 
                    % (use therefore the full MoM submatrix)
                    useACAtmp = Const.useACA;
                    Const.useACA = 0;
                    [Zact, Uact, Vact] = calcZmn(Const,zMatrices,freq,m,m,...
                        domain_m_basis_functions,domain_m_basis_functions);
                    Const.useACA = useACAtmp;

                    % NGF-en. DGFM correction here:
                    if ((Const.useDGFMmethod == 1) && (Const.runNGFenDGFMsolver))
                        % Correct the self-term here with coupling from the static domain.
                        % The other terms in the active impedance matrix eq. will be corrected individually
                        % before weighting them with the alpha-factors.
                        Zact  = Zact - Zms*Zsd_glob(:,(m-1)*Ndgfm+1:m*Ndgfm);
                        %Zact  = Zact - Zcorr((m-1)*Ndgfm+1:m*Ndgfm,(m-1)*Ndgfm+1:m*Ndgfm);
                    end%if
                    if (~Const.no_mutual_coupling_array)
                        % If we consider mutual coupling, then we add the coupling
                        % matrices to the active impedance matrix of the domain
                        % Initialisations for the ACA
                        Uact = [];
                        Vact = [];
                        for n=1:numArrayEls % N_LOOP
                            % Do not evaluate self-terms here.
                            if (m == n)
                                continue;
                            end%if

                            % Extract the coupling matrix, Zmn, scaled with alpha_nm.
                            % See issue FEKDDM-6.2: Improve the array numbering by using a
                            % bottom and top basis function index (particularly for use with
                            % the NGF-en. DGFM)
                            %domain_bot_n = Const.arrayElBasisFunctRange(n,1);
                            %domain_top_n = Const.arrayElBasisFunctRange(n,2);
                            domain_n_basis_functions = Solver_setup.rwg_basis_functions_domains{n};

                            % See issue FEKDDM-2.1: The coupling matrix can now
                            % be extracted using the ACA (call a generic getZmn routine)
                            %Zmn = zMatrices.values(domain_bot_m:domain_top_m,domain_bot_n:domain_top_n);
                            [Zmn, Umn, Vmn] = calcZmn(Const,zMatrices,freq,m,n,domain_m_basis_functions,...
                                domain_n_basis_functions);
                            if ((Const.DGFMweightVectorCalcScheme == 0) || (Const.DGFMweightVectorCalcScheme == 1) ...
                                || (Const.DGFMweightVectorCalcScheme == 5))
                                % alpha_nm is a scalar
                                alphanm = dgfm.weightVectors(n,index)/dgfm.weightVectors(m,index);                            
                            else
                                % alpha_nm is a matrix (see FEKDDM-1.2 and also notes made on 21.06.2013)
                                alphanm = complex(zeros(Ndgfm,Ndgfm));
                                for jj = 1:Ndgfm
                                    alphanm(jj,:) = (dgfm.weightVectors(:,n,index)./dgfm.weightVectors(:,m,index)).';
                                end%for
                            end%if
                            if ((Const.useDGFMmethod == 1) && (Const.runNGFenDGFMsolver))
                                % -- NGF-en. DGFM - add correction terms
                                % Correct the coupling terms here with coupling from the static domain.
                                %Zact  = Zact  + alphanm .* (Zmn - Zcorr((m-1)*Ndgfm+1:m*Ndgfm,(n-1)*Ndgfm+1:n*Ndgfm));
                                Zact  = Zact  + alphanm .* (Zmn - Zms*Zsd_glob(:,(n-1)*Ndgfm+1:n*Ndgfm));
                            else
                                % -- Standard DGFM - no correction terms
                                % ===> TO-DO: Check here that the DGFM weighting vector is calculated correctly for a frequency
                                %      loop
                                if (Const.useACA && Const.ACAalg == 3)
                                    % Concatenate the U and V rectangular matrices (see Workbook, pp. 142 - 145).
                                    if ((Const.DGFMweightVectorCalcScheme == 0) || (Const.DGFMweightVectorCalcScheme == 1))
                                        % alpha_nm is a scalar
                                        alphanm = dgfm.weightVectors(n,index)/dgfm.weightVectors(m,index);
                                    else
                                        % alpha_nm is a rectangular matrix (see FEKDDM-2.1 and also notes made on 01.04.2013, 
                                        % Workbook, pp. 143)
                                        r = size(Vmn,1); % Extract here the rank of V
                                        alphanm = complex(zeros(r,Ndgfm));
                                        for jj = 1:r
                                            alphanm(jj,:) = (dgfm.weightVectors(:,n,index)./dgfm.weightVectors(:,m,index)).';
                                        end%for
                                    end
                                    Uact = [Uact Umn];
                                    % Correct Vmn with the alpha scaling coefficients
                                    Vmn = alphanm.*Vmn;
                                    Vact = [Vact; Vmn];
                                else
                                    Zact = Zact + alphanm .* Zmn;
                                end
                                % TO-DO: Plot here the convergence of the
                                % diagonal and off diagonal terms (if we are
                                % considering a problem that is sorted
                                % according to mutual coupling, see [1])
                            end%if

                        end % N_LOOP
                    end%if(~no_coupling)

                    if (Const.useDGFMmethod == 1)
                        % -- DGFM method 1: Solve on active element level

                        if (Const.useACA && Const.ACAalg == 3)
                        % See FEKDDM-2.1: If we are using the ACA with U and V
                        % only, then we need to calculate Zact = U*V
                            Zact = Zact + Uact*Vact;
                        end

                        % Solve the active / scan impedance matrix equation for the element
                        % and store the current in the correct "global" MoM vector
                        % NOTE: Distinguish between two methods for the NGF-DGFM below
                        if (Const.runNGFenDGFMsolver)
                        % See issue FEKDDM-6.1: For the NGF-enhanced DGFM we add also a
                        % coupling here from the excitation that is defined in the static domain here
                            b = ngf.L(:,:,freq)\yVectors.values(1:Nngf,index);
                            Vstat = Zms*(ngf.U(:,:,freq)\b);
                            Varr = yVectors.values(domain_m_basis_functions,index) - Vstat;
                        else
                        % Only array environment - no coupling from static domain
                            Varr = yVectors.values(domain_m_basis_functions,index);
                        end
                        dgfm.Isol(domain_m_basis_functions,index) = Zact \ Varr;
                        
                        % Check first here whether we are using Interpolation. If so, then we need 
                        % to calculate the MBF coefficients also for this sample on the interpolation grid.
                        if (Const.useDGFMinterpolation)                        
                            % Now we can find the MBF weights by setting up a reduced impedance matrix
                            Zred = (dgfm.refMBFs)' * Zact * dgfm.refMBFs;
                            Vrwg = yVectors.values(domain_m_basis_functions,1);
                            Vred = (dgfm.refMBFs)' * Vrwg;
                            % Store these MBF coefficients now for later processing.
                            dgfm.MBF_weights(m,:) = Zred \ Vred;

                            % Compare now this solution against that obtained from the DGFM for this particular sample
                            if (true)
                                % First build again Isol using the MBFs
                                Isol_m = dgfm.refMBFs * (dgfm.MBF_weights(m,:)).';
                                relError = calculateErrorNormPercentage(dgfm.Isol(domain_m_basis_functions,1), Isol_m);
                                message_fc(Const,sprintf('Rel. error norm. for element. %d compared to DGFM sol. %f percent',m, relError));
                            end
                        end %if (Const.useDGFMinterpolation)

                    elseif(Const.useDGFMmethod == 2)
                        % -- DGFM method 2: Solve on block matrix level (after M_LOOP)
                        dgfm.Zact(domain_m_basis_functions,domain_m_basis_functions,index) = Zact;
                    end%if

                    % For improvements of the DGFM-enhanced Jacobi Method, or the iterative DGFM method
                    % we need to store also the active impedance matrices on the diagonal
                    if (Const.storeZact)
                        dgfm.Zact(domain_m_basis_functions,domain_m_basis_functions,index) = Zact;
                    end%if

                else % -- Array element m not active.
                     % Passive elements are treated as follows, based on the value of Const.calcDGFMweightVectors.m:
                     switch Const.DGFMweightVectorCalcScheme

                        case 0 % Unit weighting (1 + 1i*0), set dgfm.Isol(arrEl_m) = {0}
    			            dgfm.Isol((m-1)*Ndgfm+1:m*Ndgfm,index) = complex(zeros(Ndgfm,1));

                        case 1 % Use excitation ratio, set dgfm.Isol(arrEl_m) = {0}
    			            dgfm.Isol((m-1)*Ndgfm+1:m*Ndgfm,index) = complex(zeros(Ndgfm,1));

                        case 2 % Use the FEKO MoM Xsol (exact) solution
    			            %dgfm.Isol((m-1)*Ndgfm+1:m*Ndgfm,index) = xVectors.values(domain_bot_m:domain_top_m,index);
                            dgfm.Isol(domain_m_basis_functions,index) = xVectors.values(domain_m_basis_functions,index);
                            
                        case 3 % Use the secondary MBFs induced by the other active elements in the array,
                               % i.e. the Jacobi solution for this element as returned from the routine
                               % calcDGFMweightVectors.m (see above call)
                            dgfm.Isol((m-1)*Ndgfm+1:m*Ndgfm,index) = jackDGFM.Isol((m-1)*Ndgfm+1:m*Ndgfm,index);

                        otherwise %Internal Error - should have been checked earlier
                            error([sprintf('[runDGFMsolver] Internal Error')]);
                            mesage(Const,sprintf('[runDGFMsolver] Internal Error'));

    				end%switch
                end % if (Const.isArrayElementActive())
            end % M_LOOP

            % If we use the DGFM sol. method 2, then solve the whole block matrix
            % equation with the active impedance matrices on the diagonal.
            % See issue FEKDDM-10 and also checkConfParamsForVersion.m: The
            % following case is only supported for active arrays, i.e. where
            % all the elements are excited in the solution configuration. Check
            % already done in checkConfParamsForVersion.m.
            if (Const.useDGFMmethod == 2)
                % If we use the NGF-en. DGFM, then apply the same correction as for
                % the above case (i.e DGFM sol. method 1). This corresponds to
                % NGF-en. formulation as written down in workbook, pp. 108-109.
                if (Const.runNGFenDGFMsolver)
                    % See issue FEKDDM-6.1: For the NGF-enhanced DGFM, we need to add the correction/coupling
                    % term from the static interaction matrix here
                        Zds = zMatrices.values(Nngf+1:Nmom,1:Nngf);
                        dgfm.Zact = dgfm.Zact - Zds*Zsd_glob;
                    % See issue FEKDDM-6.1: For the NGF-enhanced DGFM we add also a
                    % coupling here from the excitation that is defined in the static domain here
                        b = ngf.L\yVectors.values(1:Nngf,index);
                        Vstat = Zds*(ngf.U\b);
                        Varr = yVectors.values(Nngf+1:Nmom,index) - Vstat;
                    else
                    % Only array environment - no coupling from static domain
                        Varr = yVectors.values(Nngf+1:Nmom,index);
                    end
                dgfm.Isol(:,index) = dgfm.Zact(:,:,index) \ yVectors.values(Nngf+1:Nmom,index);
            end%if
        end%for solNum = solStart:solEnd

        % End timing
        dgfm.solTime(freq) = toc + dgfm.ngfdgfm_init_time(freq);

        % Set the total solution time for the DGFM (all frequency points)
        dgfm.totsolTime = dgfm.totsolTime + dgfm.solTime(freq);
        dgfm.totngfdgfm_init_time = dgfm.totngfdgfm_init_time + dgfm.ngfdgfm_init_time(freq);

        % Calculate the memory usages
        % Memory usage remains the same for each frequency iteration.
        if (Const.useDGFMmethod == 1)
            % -- Zact local matrix
            
            if ((Const.DGFMweightVectorCalcScheme == 0) || (Const.DGFMweightVectorCalcScheme == 1) ...
                || (Const.DGFMweightVectorCalcScheme == 5))
                % alpha_nm is a scalar
                %    size(Zact) + alpha_nm * size(Zpq)
                dgfm.memUsage = byteSize(complex(zeros(size(Zact,1)*2,size(Zact,2)*1)));            
            else
                % alpha_nm is a matrix (see FEKDDM-1.2 and also notes made on 21.06.2013)
                %    size(Zact) + size(alpha_nm) * size(Zpq)
                dgfm.memUsage = byteSize(complex(zeros(size(Zact,1)*3,size(Zact,2)*1)));
            end%if
            
            % Note, if the weighting coefficients are calculated as alpha^1,
            % then the memory usage of the DGFM is twice that used to store the
            % local interaction matrix.
        else
            % -- Zact  global matrix
            dgfm.memUsage = byteSize(dgfm.Zact(:,:,1));
        end
        
        % Write the DGFM solution to a ASCII str file, so that it can be read
        % again by FEKO (for plotting in POSTFEKO) - only if requested (i.e. if the filename is defined)
        if (~isempty(Const.SUNEMdgfmstrfilename))
            writeSolToFile(Const, dgfm);
        end%if

    % --------------------- END HERE FREQUENCY LOOP    
    end%for freq=1:numFreq
    
    message_fc(Const,sprintf('Finished DGFM solver in %f sec.',dgfm.totsolTime));
    if (Const.runNGFenDGFMsolver)
        message_fc(Const,sprintf('NGF-en. DGFM time (Jd): %f sec.',dgfm.totngfdgfm_init_time));
    end%if
    message_fc(Const,sprintf('Memory usage of DGFM (storing Zact) %s',dgfm.memUsage));
    % TO-DO: Also add memory usage for Zmn and also for the 2 x weighting vectors        

    for freq = 1:numFreq
        for solNum=1:numRHSperFreq
        index = solNum + (freq-1)*numRHSperFreq;
            % Compare the DGFM solution obtained with MATLAB, with that obtained by FEKO
            % that was stored in xVectors.values
            % See issue FEKDDM-6.2: If the NGF-enhanced DGFM is used, then we compare
            % only the current that is on the dynamic domain, i.e. the finite array.
            % Update also the above now to account for the IFB-DGFM (general offset Const.domAoffset is
            % calculated in the function extractBasisFunctionSetup.m). When we bring back NGF-DGFM, then check from where the error is calculated.
            %dgfm.relError(index) = calculateErrorNormPercentage(xVectors.values(Const.domAoffset+1:Nmom,index), dgfm.Isol(:,index));
            dgfm.relError(index) = calculateErrorNormPercentage(xVectors.Isol(1:Nmom,index), dgfm.Isol(:,index));
            message_fc(Const,sprintf('Rel. error norm. for Sol. %d of %d of freq. %d of %d compared to FEKO sol. %f percent',solNum, numSols, ...
                freq, numFreq, dgfm.relError(index)));
        end
    end%for freq = 1:numFreq


