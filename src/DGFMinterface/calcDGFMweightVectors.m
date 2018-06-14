function [weightVectors, jackDGFM] = calcDGFMweightVectors(Const, Solver_setup, zMatrices, yVectors, xVectors, mbfs)
    %calcDGFMweightVectors v0.2
    %   Date: 28.11.2013
    %   Usage:
    %       [weightVectors] = calcDGFMweightVectors(Const, zMatrices, yVectors, xVectors)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
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
    %       weightVectors
    %           An (numBasisf x numArrayEls) matrix, where each collumn is
    %           the associated complex coefficients that is used to
    %           calculate the alphamn weighting vector when the active /
    %           scan impedance matrix of the DGFM is assembled.
    %       jackDGFM
    %           The Jacobi iterations as called here (i.e. internal to
    %           the DGFM, for the treatment of passive arrays (see
    %           FEKDDM-2.5)

    %   Description:
    %       Determines the complex coefficients that is used to
    %       calculate the alphamn weighting vector when the active /
    %       scan impedance matrix of the DGFM is assembled.
    %       See FEKDDM-5: Added now multiple solution configurations,
    %       with the benefit of treating passive array environments, e.g.
    %       for the calculation of embedded element patterns.
    %
    %   =======================
    %   Written by Danie Ludick on June 21, 2013.
    %   Last updated on NOvember 28, 2013.
    %   EMSS-SA (Pty) Ltd
    %   Email: dludick@emss.co.za

    %   References
    % [1] "Efficient Analysis of Large Irregular Antenna Arrays using the
    %      Domain Green's Function Method", Special Issue on Antennas &
    %      Propagation Submission, 2013, D.J. Ludick, R. Maaskant, et. al

    narginchk(6,6);

    message_fc(Const,sprintf('  Calculating the DGFM weighting vectors'));

    Nmom = Solver_setup.num_mom_basis_functions;   % Total number of basis functions for whole problem
    %Nngf = Solver_setup.num_ngf_basis_functions;   % Number of basis functions for NGF domain
    %NtotArray = Nmom;                              % Total number of basis functions for the array
    Ndgfm = Solver_setup.mom_basis_functions_per_array_element;  % Number of basis functions per array element
    numArrayEls = Solver_setup.num_finite_array_elements;        % The number of array elements
    numSols = 1;                                   % The number of solutions configurations (TO-DO: only have 1 for now)

    % We need to run the DGFM (or i-DGFM) for each frequency point.
    numFreq = Solver_setup.frequencies.freq_num;   % The number of frequency points to process
    numRHSperFreq = numSols / numFreq;            % The number of solutions per frequency point

    % Some info about the solution configurations. Already printed out in
    % runDGFMsolver.m
    %message_fc(Const,sprintf('  numSols : %d', numSols));
    %message_fc(Const,sprintf('  numFreq : %d', numFreq));
    %message_fc(Const,sprintf('  numRHSperFreq : %d', numRHSperFreq));

    % Initialise here the return values before we access the frequency loop
    if ((Const.DGFMweightVectorCalcScheme == 0)||(Const.DGFMweightVectorCalcScheme == 1)||...
        (Const.DGFMweightVectorCalcScheme == 5))
        weightVectors  = complex(ones(numArrayEls,numSols));
    elseif ((Const.DGFMweightVectorCalcScheme == 2) || (Const.DGFMweightVectorCalcScheme == 3)||...
            (Const.DGFMweightVectorCalcScheme == 4))
        weightVectors  = complex(ones(Ndgfm,numArrayEls,numSols));
    end%if

    % The weighting vectors are associated with the different frequency iterations.
    for freq=1:numFreq

        % Extract the solutions from which to start and end our CBFM process
        solStart = 1;
        solEnd   = numRHSperFreq;

        % Initialise the return values
        %weightVectors  = [];
        jackDGFM = [];

        switch Const.DGFMweightVectorCalcScheme

            case 0 % Unity - a scalar multiplication done in runDGFMsolver.m
                
                for solNum = solStart:solEnd
                    % Note: for frequency loops we need to access the correct element
                    index = solNum + (freq-1)*numRHSperFreq;
                    for arrEl = 1:numArrayEls
                        % Only use unity weight when the element is active,
                        % otherwise it is inactive and it gets a weight of zero
                        if (Const.is_array_element_active(arrEl,index))
                            weightVectors(arrEl,index) = 1 + 1i*0;
                        else
                            weightVectors(arrEl,index) = 0 + 1i*0;
                        end
                    end%for
                end%for

            case 1 % Standard - use the ratio of the applied excitations
                
                % Loop over all the array elements and copy the corresponding excitation values
                % for the elements. NOTE: Currently we allow only localized excitations, i.e. where
                % e.g. edge feeds for strip dipoles - or no sources at all, i.e. passive array elements.
                for solNum = solStart:solEnd
                    % Note: for frequency loops we need to access the correct element
                    index = solNum + (freq-1)*numRHSperFreq;
                    for arrEl = 1:numArrayEls
                        domain_bot = Const.arrayElBasisFunctRange(arrEl,1); % First el. bottom
                        domain_top = Const.arrayElBasisFunctRange(arrEl,2); % Last el. top
                        
                        sourcePos_arrEl  = find(yVectors.values(domain_bot:domain_top,index));
                        numSources_arrEl = length(sourcePos_arrEl);

                        if ((numSources_arrEl ~= 1) && (numSources_arrEl ~= 0) && ((numSources_arrEl ~= Ndgfm)))
                            message_fc(Const,'[calcDGFMweightVectors] Invalid source specification for DGFM Weighting Scheme 1');
                            error ('[calcDGFMweightVectors] Invalid source specification for DGFM Weighting Scheme 1');
                        elseif (numSources_arrEl == Ndgfm)
                        % Add now support for distributed / global sources,
                        % such as plane waves, herzian dipoles, etc.
                            message_fc(Const,'[calcDGFMweightVectors] Extracting from plane wave source not yet supported');
                            error ('[calcDGFMweightVectors] Extracting from plane wave source not yet supported');
                            % TO-DO: Danie, the following is not working correctly (need to extract unit magnitude and
                            % correct phase ratio based on the element position). The following is not working:
                            weightVectors(arrEl,index) = sum(yVectors.values(sourcePos_arrEl+domain_bot-1,index));
                        end%if

                        if (numSources_arrEl == 0)
                            % Essentially a passive element - i.e. with a zero weighting
                            weightVectors(arrEl,index) = 0 + 1i*0;
                        elseif (numSources_arrEl == 0)
                            % Local sources (distributed sources, e.g. plane
                            % waves, already set above)
                            weightVectors(arrEl,index) = yVectors.values(sourcePos_arrEl+domain_bot-1,index);
                        end
                    end%for
                end%for


            case 2 % Use the FEKO MoM Xsol (exact) solution
                % Loop over all the array elements and copy the corresponding Xsol values
                % See FEKDDM-10: Added now support for multiple solution
                % configurations, each essentially resulting in a weighting
                % vector for each of the array elements.
                for solNum = solStart:solEnd
                    % Note: for frequency loops we need to access the correct element
                    index = solNum + (freq-1)*numRHSperFreq;
                    for arrEl = 1:numArrayEls
                        %domain_bot = Const.arrayElBasisFunctRange(arrEl,1); % First el. bottom
                        %domain_top = Const.arrayElBasisFunctRange(arrEl,2); % Last el. top
                        domain_basis_functions = Solver_setup.rwg_basis_functions_domains{arrEl};                        
                        weightVectors(:,arrEl,index) = (xVectors.Isol(domain_basis_functions,index)).';
                    end%for
                end%for

            case 3 % Use the Jacobi-Generated CBFs (see [1], Appendix A)

                if (numFreq ~= 1)
                    message_fc(Const,'[calcDGFMweightVectors] Only 1 freq. allowed for this DGFM weighting scheme');
                    error ('[calcDGFMweightVectors] Only 1 freq. allowed for this DGFM weighting scheme');
                end%if

                % Set a constant to suppress / format some output in the JACKIT solver
                Const.runJACKITfromDGFM = true;
                [jackDGFM] = runJACKITsolver(Const, Solver_setup, zMatrices, yVectors, xVectors, mbfs);
                Const.runJACKITfromDGFM = false;
                % Loop over all the array elements and copy the corresponding Xsol values
                % See FEKDDM-10 (and above comment): Added now support for multiple solution
                % configurations, each essentially resulting in a weighting
                % vector for each of the array elements.
                % TO-DO: The above has to be expanded for multiple frequency loops. The loop below takes
                % this into account, but we need to add support for numFreq>1 inside runJACKITsolver()
                for solNum = solStart:solEnd
                    % Note: for frequency loops we need to access the correct element (see above comment)
                    index = solNum + (freq-1)*numRHSperFreq;
                    for arrEl = 1:numArrayEls
                        weightVectors(:,arrEl,index) = jackDGFM.Isol((arrEl-1)*Ndgfm+1:arrEl*Ndgfm,index);
                    end%for
                end%for

            case 4 % Use an arbitrary excitation vector (useful for use with distributed sources, or
                   % for IFB related runs where the RHS changes rapidly)
                % Loop over all the array elements and copy the corresponding yVector values. We need
                % to account also for the fact that the values may be offset (use therefore the )
                for solNum = solStart:solEnd
                    index = solNum + (freq-1)*numRHSperFreq;
                    for arrEl = 1:numArrayEls
                        % See issue FEKDDM-6.2: Improve the array numbering by using a
                        % bottom and top basis function index (particularly for use with
                        % the NGF-en. DGFM) - also for the IFB-DGFM (where the finite array is offset in the
                        % matrix)
                        domain_bot = Const.arrayElBasisFunctRange(arrEl,1);
                        domain_top = Const.arrayElBasisFunctRange(arrEl,2);
                        % Update the weight vectors (need to use the i-DGFM)
                        weightVectors(:,arrEl,index) = yVectors.values(domain_bot:domain_top,index);
                    end%for
                end%for

            case 5 % Calculate the sum of the currents on the array elements (stored in xVector.values)
                   % - a scalar multiplication done in runDGFMsolver.m. This is used in the iterative-
                   % DGFM calculation
                for solNum = solStart:solEnd
                    index = solNum + (freq-1)*numRHSperFreq;
                    for arrEl = 1:numArrayEls
                        % Use the sum or the average of the currents on the array element
                        domain_bot = Const.arrayElBasisFunctRange(arrEl,1); % First el. bottom
                        domain_top = Const.arrayElBasisFunctRange(arrEl,2); % Last el. top
                        weightVectors(arrEl,index) = mean(xVectors.values(domain_bot:domain_top,index));
                    end%for
                end%for

            otherwise %Error
                error([sprintf('Case %d not supported for calculating DGFM weight vectors', ...
                        Const.DGFMweightVectorCalcScheme)]);
                mesage(Const,sprintf('Case %d not supported for calculating DGFM weight vectors',...
                    Const.DGFMweightVectorCalcScheme));

    	end%switch

    end%for freq=1:numFreq

    message_fc(Const,sprintf('  Finished calculating the DGFM weighting vectors'));