function [cbfm] = runCBFMsolver(Const, Solver_setup, zMatrices, yVectors, xVectors, mbfs)
    %runCBFMsolver
    %   Date: 18.08.2015
    %   Usage:
    %       [cbfm] = runCBFMsolver(Const, Solver_setup, zMatrices, yVectors, xVectors, mbfs)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details    
    %       zMatrices
    %           The Z-matrices data
    %       yVectors
    %           The Yrhs-vector data
    %       xVectors
    %           The Xsol-vector data (i.e. MoM solution of FEKO)
    %       mbfs
    %           The MBFs (prim & sec.) generated in runMBFgenerator
    %
    %   Output Arguments:
    %       cbfm
    %           Structs containing CBFM solution and timing data
    %
    %   Description:
    %       [~2013]
    %         Runs the Characteristic Basis Function Method (CBFM) solution based
    %         on the Z and Y data that was read / parsed from the FEKO *.out,
    %         *.mat, *.str and *.rhs filesand also the array data that was setup
    %       [~2015]
    %         Also added now the use of the reduced and orthonormalised MBFs (as done via
    %         the SVD approach in the MBF generator routine).
    %
    %   References:
    %   [1] V. V. S. Prakash and Raj Mittra, "Characteristic Basis Function Method:
    %       A New Technique for Efficient Solution of Method of Moments Matrix Equations,"
    %       in Microwave and Optical Technology Letters, Vol. 36, No. 2
    %       Jan, 2003, pp. 95-100.
    %   [2] R. Maaskant, D. Ludick, C. Bencivenni, M. V. Ivashina, D. Davidson,
    %       "A Compressed Sensing Technique Applied to the Characteristic Basis Function Method"
    %       2015 [to be published]
    %   =======================
    %   Written by Danie Ludick on June 21, 2013 and updated in 2015 and 2018, respectively.
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    narginchk(6,6);

    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running CBFM solver'));
    if (Const.no_mutual_coupling_array)
        message_fc(Const,sprintf('(*** Mutual coupling between array elements cannot be ignored ***)'));
        message_fc(Const,sprintf('(*** for the CBFM, i.e. primary CBFs will still be used to    ***)'));
        message_fc(Const,sprintf('(*** calculate the reduced impedance matrix equation          ***)'));
    end%if
    if (Const.useACA)
        message_fc(Const,sprintf('(Using the ACA)'));
        if (Const.ACAalg == 3)
            message_fc(Const,sprintf('(efficient ACA algorithm - UV only)'));
        elseif (Const.ACAalg == 2)
            message_fc(Const,sprintf('(efficient ACA algorithm)'));
        else
            message_fc(Const,sprintf('(non-efficient ACA algorithm)'));
        end
    end%if

    if (Const.useMBFreduction)
        message_fc(Const,sprintf('(Using the (SVD) reduced MBF set)'));
    end%if

    % Initialisations
    cbfm  = [];
    cbfm.name = 'cbfm';
    Nmom = Solver_setup.num_mom_basis_functions;                   % Total number of basis functions for whole problem
    Nngf = Solver_setup.num_ngf_basis_functions;                   % Number of basis functions for NGF domain
    %Ndom = Solver_setup.max_mom_basis_functions_per_array_element; % Total number of basis functions for the array
    numArrayEls = Solver_setup.num_finite_array_elements;          % The number of array elements
    
    numSols = xVectors.numSols;                % The number of reference solutions
    cbfm.numSols = numSols;                    % Calculate a solution for each configuration
    cbfm.Isol = complex(zeros(Nmom,numSols));  % The solved expansion coefficients (global)

    % Flag to control whether the relative residuum will be calculated or not (computationally expensive)
    calculateRelativeResiduum = false;
    
    % Note: We assume here that the MBFs have already been generated prior
    % to this part of the CBFM (i.e. the reduce matrix setup + solution)
    if (isempty(mbfs))
        message_fc(Const,'[runCBFMsolver] No MBFs have been setup');
        error ('[runCBFMsolver] No MBFs have been setup');
    end%if

    % Extract the solutions from which to start and end our CBFM process
    solStart = Const.solStart;
    solEnd   = Const.solEnd;

    % See FEKDDM-10: For the CBFM we need to follow the standard
    % formulation as discussed in [1], i.e. where the Zred is calculated
    % based on the active array, i.e. using all the MBFs (both primary and
    % secondary), as setup during the MBF generation phase for the case
    % where all the array elements are excited. Now determine for which
    % excitation law, all the MBFs are available (typically associated with
    % the reference/full-wave solution). If no such solution is available,
    % then we need to call again the MBF-generator with an all excited
    % array.
    actSol = 0;
    for solNum = 1:numSols
        % Loop over all the elements in each solution configuration and
        % check whether they are all active
        allActive = true;
        for el = 1:numArrayEls
            if (~Const.is_array_element_active(el,solNum))
                % Passive element detected
                allActive = false;
            end
        end%for
        if (allActive)
            % We found a solution
            actSol = solNum;
        end%if
    end

    % TO-DO: If no active solution could be found, then we need to recalculate the
    % MBFs on the array for an all-active case. For now, this is an error.
    if (actSol == 0)
        message_fc(Const,'[runCBFMsolver] For the CBFM, we need at least one active array solution configuration');
        error ('[runCBFMsolver] For the CBFM, we need at least one active array solution configuration');
    end

    % Start timing
    tic

    % See issue FEKDDM-10: We added now support for multiple solution
    % configurations. Each of the solution configurations has given rise to
    % a suitable weightVectors (based on the excitation-law, exact current,
    % etc. for that specific solution). We need to repeat the DGFM active
    % impedance matrix calculation of each array element for each solution
    % configuration.
    % See also FEKDDM-3.9:
    Zred_not_calculated = true;
    allocated = false;
    for solNum = solStart:solEnd

        if (Zred_not_calculated)
            % Allocate space for the CBFM reduced impedance matrix equation (Zred and Vred)
            % This is dependent on the different approaches for which the MBFs have been generated
            % (i.e. by using the SVD or not) or the solution we apply, i.e. CBFM vs. CS-CBFM.
            if (Const.useMBFreduction)
                % Each domain can have its own number of (reduced) primaries / secondaries. This
                % was already calculated in the MBF generator when the reduced MBFs were setup
                numRows = mbfs.totRedMBFs;
                numCols = mbfs.totRedMBFs;
            elseif (Const.useCSCBFM)
                % See [2] : Create an undetermined system (more collumns than rows)
                % Use the primary and secondary MBFs (i.e. not the SVD reduced subset) for the basis
                % functions (columns)
                numCBFsperDomain = mbfs.numPrimMBFs(1,actSol) + mbfs.numSecMBFs(1,actSol);
                numCols = numArrayEls*numCBFsperDomain;
                % Use only primary MBFS (i.e. not the SVD reduced subset)
                % for the testing functions (rows)
                numCBFsperDomain = mbfs.numPrimMBFs(1,actSol);
                numRows = numArrayEls*numCBFsperDomain;
            else % Conventional CBFM as in [1] assuming primaries + secondaries
                % See [1] : Assume here the same number of primaries and seconfaries / domain
                numCBFsperDomain = mbfs.numPrimMBFs(1,actSol) + mbfs.numSecMBFs(1,actSol);
                numRows = numArrayEls*mbfs.numPrimMBFs(1,actSol);
                numCols = numArrayEls*numCBFsperDomain;
            end

            % Note, only allocate Zred and Vred once
            if (~allocated)
                cbfm.Zred = complex(zeros(numRows, numCols));
                cbfm.Vred = complex(zeros(numRows, 1)); % Will be recalculated for each solution configuration
                allocated = true;
            end

            % ======================================
            % Setup the reduced impedance matrix Zred (done only once)
            % ======================================

            % -- Domain P
            for p=1:numArrayEls % Loop over rows
                
               % Extract basis function indices associated with domains p
               domain_p_basis_functions = Solver_setup.rwg_basis_functions_domains{p};
               Ndom_p = length(domain_p_basis_functions);

                % First create a collumn augmented matrix containing all the MBFs of domain p,
                % i.e. both primary and secondary MBFs:
                if (Const.useMBFreduction)
                    % Copy only the reduced MBFs
                    domainPcbfs = complex(zeros(Ndom_p,mbfs.numRedMBFs(p,actSol)));
                    for red_ii=1:mbfs.numRedMBFs(p,actSol)
                        domainPcbfs(:,red_ii) = mbfs.RedIsol(domain_p_basis_functions,red_ii,p,actSol);
                    end%for mbfs.numRedMBFs(p)

                elseif (Const.useCSCBFM)
                    % Now as explained in [2], for the testing MBFs we use only the primaries
                    domainPcbfs = complex(zeros(Ndom_p,mbfs.numPrimMBFs(p,actSol)));
                    for prim_ii=1:mbfs.numPrimMBFs(p,actSol)
                        domainPcbfs(:,prim_ii) = mbfs.PrimIsol(domain_p_basis_functions,prim_ii,p,actSol);
                    end%for mbfs.numPrimMBFs(p)
                else
                    % Conventional method as explained in [1] - primaries and secondaries
                    domainPcbfs = complex(zeros(Ndom_p,mbfs.numPrimMBFs(p,actSol) + mbfs.numSecMBFs(p,actSol)));
                    for prim_ii=1:mbfs.numPrimMBFs(p,actSol)
                        domainPcbfs(:,prim_ii) = mbfs.PrimIsol(domain_p_basis_functions,prim_ii,p,actSol);
                    end%for mbfs.numPrimMBFs(p)
                    for sec_ii=1:mbfs.numSecMBFs(p,actSol)
                        domainPcbfs(:,mbfs.numPrimMBFs(p,actSol) + sec_ii) = mbfs.SecIsol(domain_p_basis_functions,sec_ii,p,actSol);
                    end%for mbfs.numPrimMBFs(p)
                end

                % -- Domain Q
                for q=1:numArrayEls % Loop over collumns
                                        
                    % Do the same as the above here for domain Q:
                    
                    % Extract basis function indices associated with domains q
                    domain_q_basis_functions = Solver_setup.rwg_basis_functions_domains{q};
                    Ndom_q = length(domain_q_basis_functions);
                    
                    if (Const.useMBFreduction)
                        % Copy only the reduced MBFs
                        domainQcbfs = complex(zeros(Ndom_q,mbfs.numRedMBFs(q,actSol)));
                        for red_ii=1:mbfs.numRedMBFs(q,actSol)
                            domainQcbfs(:,red_ii) = mbfs.RedIsol(domain_q_basis_functions,red_ii,q,actSol);
                        end%for mbfs.numRedMBFs(p)

                    else % Conventional + CS-CBFM (use both primaries and secondaries)
                        % First create a collumn augmented matrix containing all the MBFs of domain q,
                        % i.e. both primary and secondary MBFs:
                        domainQcbfs = complex(zeros(Ndom_q,mbfs.numPrimMBFs(q,actSol) + mbfs.numSecMBFs(q,actSol)));
                        for prim_ii=1:mbfs.numPrimMBFs(q,actSol)
                            domainQcbfs(:,prim_ii) = mbfs.PrimIsol(domain_q_basis_functions,prim_ii,q,actSol);
                        end%for
                        for sec_ii=1:mbfs.numSecMBFs(q,actSol)
                            domainQcbfs(:,mbfs.numPrimMBFs(q,actSol) + sec_ii) = mbfs.SecIsol(domain_q_basis_functions,sec_ii,q,actSol);
                        end%for
                    end
                                        
                    % Calculate now Zred = (Jp)^t * Zpq * Jq (where t is the transpose operator)
                    % See issue FEKDDM-3.2: Add now ACA support for the coupling matrix assebly here (only for off-diagonal terms)
                    freq = 1;
                    if (p~=q)
                        [Zcoupl, Ucoupl, Vcoupl] = calcZmn(Const,zMatrices,freq,p,q,domain_p_basis_functions,domain_q_basis_functions);
                    else
                        useACAtmp = Const.useACA;
                        Const.useACA = 0;
                        [Zcoupl, Ucoupl, Vcoupl] = calcZmn(Const,zMatrices,freq,p,q,domain_p_basis_functions,domain_q_basis_functions);
                        %Zcoupl = zMatrices.values((p-1)*Ndom+1:p*Ndom,(q-1)*Ndom+1:q*Ndom);
                        Const.useACA = useACAtmp;
                    end

                    % ------------------------------------------------------------------------------
                    % Calculate the offsets for storing the Zred submatrix (Zpq)

                    % ==========
                    % Range of P (rows, i.e. testing functions)
                    if (Const.useMBFreduction)
                        % Use the SVD reduced set of MBFs
                        numMBFsP = 0;
                        for domain = 1:(p-1)
                            numMBFsP = numMBFsP + mbfs.numRedMBFs(domain,actSol);
                        end%for
                        PindxStart = numMBFsP+1;
                        PindxEnd   = (PindxStart - 1) + mbfs.numRedMBFs(p,actSol);

                    elseif (Const.useCSCBFM)
                        % See [2] -  use only the primaries for domain P
                        numMBFsP = 0;
                        for domain = 1:(p-1)
                            numMBFsP = numMBFsP + mbfs.numPrimMBFs(domain,actSol);
                        end%for
                        PindxStart = numMBFsP+1;
                        PindxEnd   = (PindxStart - 1) + mbfs.numPrimMBFs(p,actSol);
                    else
                        % Conventional as in [1] - use both primaries + secondaries
                        numMBFsP = 0;
                        for domain = 1:(p-1)
                            numMBFsP = numMBFsP + mbfs.numPrimMBFs(domain,actSol) + mbfs.numSecMBFs(domain,actSol);
                        end%for
                        PindxStart = numMBFsP+1;
                        PindxEnd   = (PindxStart - 1) + mbfs.numPrimMBFs(p,actSol) + mbfs.numSecMBFs(p,actSol);
                    end%if

                    % ==========
                    % Range of Q (columns, i.e. basis functions)
                    if (Const.useMBFreduction)
                        % Use the SVD reduced set of MBFs
                        numMBFsQ = 0;
                        for domain = 1:(q-1)
                            numMBFsQ = numMBFsQ + mbfs.numRedMBFs(domain,actSol);
                        end%for
                        QindxStart = numMBFsQ+1;
                        QindxEnd   = (QindxStart - 1) + mbfs.numRedMBFs(q,actSol);

                    else
                        % Conventional as in [1] or CS-CBFM as in [2] - use both primaries + secondaries
                        numMBFsQ = 0;
                        for domain = 1:(q-1)
                            numMBFsQ = numMBFsQ + mbfs.numPrimMBFs(domain,actSol) + mbfs.numSecMBFs(domain,actSol);
                        end%for
                        QindxStart = numMBFsQ+1;
                        QindxEnd   = (QindxStart - 1) + mbfs.numPrimMBFs(q,actSol) + mbfs.numSecMBFs(q,actSol);
                    end%if

                    % ------------------------------------------------------------------------------
                    % Block assignment for submatrix in Zred.
                    if (~Const.useACA)
                        % -- Standard (no ACA)
                        cbfm.Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = (domainPcbfs)' * Zcoupl * domainQcbfs;
                    else
                        % -- Use ACA - we have U and V matrices for off-diagonal entries
                        if (p~=q)
                            % -- Use factored U and V matrices to speed-up the matrix-vector products
                            cbfm.Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = (domainPcbfs)' * Ucoupl * Vcoupl * domainQcbfs;  
                        else
                            % Diagonal entries calculated using standard approach (Zpq)
                            cbfm.Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = (domainPcbfs)' * Zcoupl * domainQcbfs;
                        end
                    end

                end%for q=1:numArrayEls

            end%for p=1:numArrayEls

            % Debug output (NOTE: can use a lot of diskspace for many sol.
            % configurations)
            if (Const.debug)
                cd(Const.ProjectPath);
                cd([Const.OutputDirName '/Temp']);
                fid=fopen('Zred.txt','w');
                % Write out the value of Zred to a Zred.txt file that is saved in
                % the /temp directory
                for ii = 1:length(cbfm.Zred(:,1))
                    for jj = 1:length(cbfm.Zred(1,:))
                        fprintf(fid,'Zred(%d,%d) : %0.5E + j%0.5E \n',ii,jj,real(cbfm.Zred(ii,jj)),imag(cbfm.Zred(ii,jj)));
                    end%for
                end%for
                fclose(fid);

                % Also write out only the diagonal
                fid=fopen('Zred_diag.txt','w');
                for ii=1:size(cbfm.Zred(:,1))
                    fprintf(fid,'Zred_diag(%d,%d) = %0.5E + j%0.5E \n',ii,ii,real(cbfm.Zred(ii,ii)),imag(cbfm.Zred(ii,ii)));
                end%for

                fclose(fid);
                cd(Const.ProjectPath);
            end%if

            % Only calculate the LU-factorisation of the Zred matrix if we are not using the CS-CBFM
            % method (that solves the underdetermined system using convex optimisation)
            if (~Const.useCSCBFM)
                % LU-decomposition of the Z-matrix (now included here in this time)
                [L,U] = lu(cbfm.Zred);
                % Clear the Zred matrix after we have calculate the LU
                % decomposition to free some memory. Before doing so, we record
                % the size of Zred:
                cbfm.memUsage = byteSize(cbfm.Zred(:,:));
                clear('cbfm.Zred');
            end%if

            % See issue FEKDDM-3.9: For Phase array problems we recalculate
            % the reduced impedance matrix for each phase angle. (One can
            % also perhaps alternatively use a plane wave spectrum, but
            % this adds to the complexity of Zred filling as well - then
            % more primary CBFs)
            if (Const.isPhasedArray)
                Zred_not_calculated = true;
                allocated = false;
            end
        end%if (Zred_not_calculated)
                
        % ======================================
        % Setup the reduced excitation vector Vred
        % (See FEKDDM-10: This is done for each solution configuration).
        % ======================================

        % -- Domain P
        for p=1:numArrayEls % Loop over rows

            % Calculate now Vred = (Jp)^t * Vrwg (where t is the transpose operator)
            domain_p_basis_functions = Solver_setup.rwg_basis_functions_domains{p};
            Ndom_p = length(domain_p_basis_functions);

            % First create a collumn augmented matrix containing all the MBFs of domain p,
            % i.e. both primary and secondary MBFs:
            if (Const.useMBFreduction)
                % Copy only the reduced MBFs
                domainPcbfs = complex(zeros(Ndom_p,mbfs.numRedMBFs(p,actSol)));
                for red_ii=1:mbfs.numRedMBFs(p,actSol)
                    domainPcbfs(:,red_ii) = mbfs.RedIsol(domain_p_basis_functions,red_ii,p,actSol);
                end%for mbfs.numRedMBFs(p)

            elseif (Const.useCSCBFM)
                % Now as explained in [2], for the testing MBFs we use only the primaries
                domainPcbfs = complex(zeros(Ndom_p,mbfs.numPrimMBFs(p,actSol)));
                for prim_ii=1:mbfs.numPrimMBFs(p,actSol)
                    domainPcbfs(:,prim_ii) = mbfs.PrimIsol(domain_p_basis_functions,prim_ii,p,actSol);
                end%for mbfs.numPrimMBFs(p)
            else
                % Conventional method as explained in [1] - primaries and secondaries
                domainPcbfs = complex(zeros(Ndom_p,mbfs.numPrimMBFs(p,actSol) + mbfs.numSecMBFs(p,actSol)));
                for prim_ii=1:mbfs.numPrimMBFs(p,actSol)
                    domainPcbfs(:,prim_ii) = mbfs.PrimIsol(domain_p_basis_functions,prim_ii,p,actSol);
                end%for mbfs.numPrimMBFs(p)
                for sec_ii=1:mbfs.numSecMBFs(p,actSol)
                    domainPcbfs(:,mbfs.numPrimMBFs(p,actSol) + sec_ii) = mbfs.SecIsol(domain_p_basis_functions,sec_ii,p,actSol);
                end%for mbfs.numPrimMBFs(p)
            end
            
            % Extract the excition vector entries associated with domain p.
            Vrwg = yVectors.values(domain_p_basis_functions,solNum);

            % ==========
            % Range of P (rows, i.e. testing functions)
            if (Const.useMBFreduction)
                % Use the SVD reduced set of MBFs
                numMBFsP = 0;
                for domain = 1:(p-1)
                    numMBFsP = numMBFsP + mbfs.numRedMBFs(domain,actSol);
                end%for
                PindxStart = numMBFsP+1;
                PindxEnd   = (PindxStart - 1) + mbfs.numRedMBFs(p,actSol);

            elseif (Const.useCSCBFM)
                % See [2] -  use only the primaries for domain P
                numMBFsP = 0;
                for domain = 1:(p-1)
                    numMBFsP = numMBFsP + mbfs.numPrimMBFs(domain,actSol);
                end%for
                PindxStart = numMBFsP+1;
                PindxEnd   = (PindxStart - 1) + mbfs.numPrimMBFs(p,actSol);
            else
                % Conventional as in [1] - use both primaries + secondaries
                numMBFsP = 0;
                for domain = 1:(p-1)
                    numMBFsP = numMBFsP + mbfs.numPrimMBFs(domain,actSol) + mbfs.numSecMBFs(domain,actSol);
                end%for
                PindxStart = numMBFsP+1;
                PindxEnd   = (PindxStart - 1) + mbfs.numPrimMBFs(p,actSol) + mbfs.numSecMBFs(p,actSol);
            end%if

            % Block assignment for subvector in Zred
            cbfm.Vred(PindxStart:PindxEnd) = (domainPcbfs)' * Vrwg;

        end%for p=1:numArrayEls

        % End timing
        cbfm.Vsetup = toc;

        % Debug output
        if (Const.debug)
            cd(Const.ProjectPath);
            cd([Const.OutputDirName '/Temp']);
            fid=fopen('Vred.txt','w');
            % Write out the value of Zred to a Zred.txt file that is saved in
            % the /temp directory
            for ii = 1:length(cbfm.Vred(:))
                fprintf(fid,'Vred(%d) : %0.5E + j%0.5E \n',ii,real(cbfm.Vred(ii)),imag(cbfm.Vred(ii)));
            end%for
            fclose(fid);
            cd(Const.ProjectPath);
        end%if

        % ======================================
        % Solve the reduced matrix equation Ired = (Zred)^(-1)Vred
        % and expand the solution for Imom
        % ======================================

        if (Const.useCSCBFM)
            % See [2], solve these equations using the CS-CBFM method
            % TO-DO: Need to figure out how to apply timing here
            % tic
            cbfm.Ired = runCVXsolver(Const, cbfm.Zred, cbfm.Vred);
            % toc
        else
            % Back-wards substitution
            b = L\cbfm.Vred;
            cbfm.Ired = U\b;
        end%if

        % Expand the Ired to Imom for each of the p array elements.
        % TO-DO: For interconnected domains, we need to average the solution
        % in the overlapping region.
        for p=1:numArrayEls

            % Calculate now Vred = (Jp)^t * Vrwg (where t is the transpose operator)
            domain_p_basis_functions = Solver_setup.rwg_basis_functions_domains{p};
            Ndom_p = length(domain_p_basis_functions);

            if (Const.useMBFreduction) % We have a reduced set of basis functions on each element

                % Calculate the total number of CBFs of the previous domain
                % so that we can correctly access the correct weighting factors
                % for domain P.
                if (p > 1)
                    offset = offset + mbfs.numRedMBFs(p-1,actSol);
                else
                    offset = 0;
                end%if
                % Add the contribution of the reduced CBFs
                for red_ii=1:mbfs.numRedMBFs(p,actSol)
                    fak = cbfm.Ired(offset + red_ii);
                    cbfm.Isol(domain_p_basis_functions,solNum) = cbfm.Isol(domain_p_basis_functions,solNum) + ...
                        fak.*mbfs.RedIsol(domain_p_basis_functions,red_ii,p,actSol);
                end%for n=1:numArrayEls

            else % use the conventional CBFM or CS-CBFM where in both cases we have primaries and secondaries

                % Calculate the total number of CBFs of the previous domain
                % so that we can correctly access the correct weighting factors
                % for domain P.
                if (p > 1)
                    offset = offset + mbfs.numPrimMBFs(p-1,actSol)+mbfs.numSecMBFs(p-1,actSol);
                else
                    offset = 0;
                end%if

                % Add the contribution of the primary CBFs
                for prim_ii=1:mbfs.numPrimMBFs(p,actSol)
                    % Note, assumed here for mbfs.PrimIsol, is that each of the
                    % domains only have 1 primary MBF. The last term in the
                    % following assignment has to be changed to: mbfs.PrimIsol(:,m,p)
                    fak = cbfm.Ired(offset + prim_ii);
                    cbfm.Isol(domain_p_basis_functions,solNum) = cbfm.Isol(domain_p_basis_functions,solNum) + ...
                        fak.*mbfs.PrimIsol(domain_p_basis_functions,prim_ii,p,actSol);
                end%for n=1:numArrayEls
                if (~Const.calcSecMBFs)
                    % If there are no Secondary MBFs, then continue with the
                    % next array element (only primary MBFs then taken into
                    % account)
                    continue;
                end%if

                %  Add the contribution of the secondary CBFs
                for sec_ii=1:mbfs.numSecMBFs(p,actSol)
                    % See FEKDDM-4.1 and also comment in runMBFgenerator: The negative sign has already been
                    % accounted for when the secondary MBFs are generated.
                    fak = cbfm.Ired(offset + mbfs.numPrimMBFs(p,actSol) + sec_ii);
                    cbfm.Isol(domain_p_basis_functions,solNum) = cbfm.Isol(domain_p_basis_functions,solNum) + ...
                        fak.*mbfs.SecIsol(domain_p_basis_functions,sec_ii,p,actSol);
                end%for
            end%if
        end%for p=1:numArrayEls        

    end%for solNum = 1:numSols

    % End timing
    cbfm.solTime = toc;

    % The CBFM memory usage is primarily dominated by the reduced impedance matrix:
    % This is now listed above, as we need to clear cbfm.Zred in order to
    % make space for the LU decomposition
    cbfm.memUsage = byteSize(cbfm.Zred(:,:));

    % Write the CBFM solution to a ASCII str file, so that it can be read
    % again by FEKO (for plotting in POSTFEKO) - only if requested (i.e. if the filename is defined)
    if (~isempty(Const.SUNEMcbfmstrfilename))
        writeSolToFile(Const, cbfm);
    end%if

    % Add also here the cost for the MBF calculation
    message_fc(Const,sprintf('Finished CBFM solver in %f sec. (including MBF calculation : %f sec.)',cbfm.solTime, ...
        cbfm.solTime + mbfs.totTime));
    % Update now the CBFM time:
    cbfm.solTime = cbfm.solTime + mbfs.totTime;
    message_fc(Const,sprintf('Memory usage of CBFM %s',cbfm.memUsage));

    for solNum = solStart:solEnd
        % Compare the CBFM solution obtained with MATLAB, with that obtained by FEKO
        % that was stored in xVectors.values
        % See issue FEKDDM-6.2: If the NGF-enhanced CBFM is used, then we compare
        % only the current that is on the dynamic domain, i.e. the finite array.
        cbfm.relError(solNum) = calculateErrorNormPercentage(xVectors.Isol(:,solNum), cbfm.Isol(:,solNum));
        message_fc(Const,sprintf('Rel. error norm. for Sol. %d of %d compared to FEKO sol. %f percent.',solNum, numSols, cbfm.relError(solNum)));

        % =======================
        % New addition: 2017-05-15: - Calculate the relative residuum
        % -------------------------
        if (calculateRelativeResiduum)

            % % Calculate the relative residual error here (only if enabled, as this is expensive)
            % if (calculateRelativeResiduum)
            %     ifbmom.relResError(k) = (calculateErrorNormPercentage(yVectors.values, zMatrices.values*ifbmom.Isol))/100.0;
            % end%if

            cbfm.relResError(solNum) = calculateErrorNormPercentage(yVectors.values(:,solNum), zMatrices.values*cbfm.Isol(:,solNum))/100.0;
            message_fc(Const,sprintf('Rel. residuum error for Sol. %d of %d compared to FEKO sol. %f',solNum, numSols, cbfm.relResError(solNum)));

            % Plot now the relative residual as a function of the runtime:
            x1 = cbfm.solTime;
            y1=cbfm.relResError;
            xlab = 'Total runtime [seconds]';
            ylab = 'Relative residuum (\epsilon)';
            % Build the title string here.
            titleString = sprintf('CBFM : useMBFreduction : %d; MBFthreshold : %f',Const.useMBFreduction, Const.MBFthreshold);
            % Build here the file-names to reflect the particular solution configuration
            img_fcd_filename = sprintf('CMA_svdthresh_%f',Const.MBFthreshold);
            imgString = strcat('relresvstime_',img_fcd_filename);
            fcdString = imgString; % Name the *.fcd name, same as the image file name.
            Const.plotSemiLogY=true;  
            plotData1(Const,x1,y1,xlab,ylab,titleString,imgString,fcdString);
            Const.plotSemiLogY=false;
        end%if

    end
