function [mbfs] = runMBFgenerator(Const, Solver_setup, zMatrices, yVectors, xVectors)
    %runMBFgenerator
    %   Date: 30.11.2013
    %   Usage:
    %       [mbfs] = runMBFgenerator(Const, zMatrices, yVectors, xVectors)
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
    %   Output Arguments:
    %       dgfm
    %           Structs containing MBFs (pimary and possibly secondaries) and timing data
    %
    %   Description:
    %       Calculates primary and secondary MBFs
    %        - Primary MBFs are generated as Jprim = (Zself)^(-1) * Vself
    %        - Secondary MBFs are generated as Jsec = (Zself)^(-1) * Zcoupl * Jprim
    %          where Zself and Zcoupl are the self-interaction and coupling matrices
    %          between two domains P and Q. Vself is the excitation vector entries
    %          local to the domain P. 
    %
    %   TO-DO: 
    %        - The Const.secMBFcouplingrad variable, specifies the coupling radius that is 
    %          used to specify which elements should be included for secondary CBFs. Currently, 
    %          all elements are included. The FEKO geometry should then be read also from the 
    %          *.out file. - logged in FEKDDM-3.4
    %        - The ACA can also be used here to accelerate the matrix vector product involved 
    %          with the secondary MBF generation - logged in FEKDDM-3.2
    %        - Parallelisation of the solver - logged in FEKDDM-3.3
    %
    %   Assumptions:
    %        - All domains are the same size (i.e. contains the same number of unknowns)
    %
    %   References:
    %   [1] V. V. S. Prakash and Raj Mittra, "Characteristic Basis Function Method: 
    %       A New Technique for Efficient Solution of Method of Moments Matrix Equations," 
    %       in Microwave and Optical Technology Letters, Vol. 36, No. 2
    %       Jan, 2003, pp. 95-100.
    %   [2] Maaskant, et. al, "Fast Analysis of Large Antenna Arrays Using the Characteristic 
    %       Basis Function Method and the Adaptive Cross Approximation Algorithm", IEEE TAP
    %
    %   =======================
    %   Written by Danie Ludick on June 24, 2013.    
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za
    %
    % Notes:
    %       See issue FEKDDM-3.5 (and FEKDDM-10): Multiple RHS vectors are now supported
    
    narginchk(5,5);
    
    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running MBF generator'));
    if (Const.no_mutual_coupling_array)
        % TO-DO: This is not correctly calculated for the "no coupling" case. We should just relay the
        % primary MBFs straight through to the answer.
        message_fc(Const,sprintf('(*** Mutual coupling between array elements ignored ***)'));
        message_fc(Const,sprintf('(*** (only accounting for primary MBFs)                ***)'));
        % Simple - just switch off calculating secondary MBFs
        Const.calcSecMBFs = false;
    end%if
        
    % Initialise the return values
    mbfs  = [];
    Nmom        = Solver_setup.num_mom_basis_functions;    % Total number of basis functions for whole problem
    Nngf        = Solver_setup.num_ngf_basis_functions;    % Number of basis functions for NGF domain

    % 2017.06.03: We need to work on the max number of basis functions (BFs) per element here (as we might have
    % interconnected) arrays that have different number of BFs / element.
    %Ndom        = Solver_setup.max_mom_basis_functions_per_array_element;   % Number of basis functions per array element        
    numArrayEls = Solver_setup.num_finite_array_elements;  % The number of array elements
    numGeneratingSubarrays = Solver_setup.generating_subarrays.number_of_domains; % Number of generating sub-arrays,
                                                                                  % for connected structures    
    numSols     = xVectors.numSols;             % The number of solutions configurations
        
    % If we are reducing and orthonormalising the MBFs by using the SVD,
    % then store these in a differents array than mbfs.PrimIsol and
    % mbfs.SecIsol
    if (Const.useMBFreduction)
        % Note, we are unsure at this stage how many reduced MBFs we will
        % retain. Just allocated enough space here (numArrayEls should do
        % the trick).
        % 2018.06.04: Increase now the size of RedIsol to the global MoM size (not per domain any more).
        mbfs.RedIsol = complex(zeros(Nmom,numArrayEls,numArrayEls,numSols));
        mbfs.numRedMBFs = zeros(numArrayEls,numSols); % Number of Red. MBFs / array element / solution config.
    end%if
        
    % See FEKDDM-3.1: Store the primary CBFs also in the same structure as
    % the secondary CBFs.
    % Structure to be followed:   (MBF{1:Ndom},PrimIndx,DomainInd)
    % NOTE: The PrimIndx is to account for the port excitation on the element, i.e. each port results 
    % in one primary MBF
    % See issue FEKDDM-3.5 (and FEKDDM-10), multiple solution
    % configurations are now supported (i.e. for multiple RHS vectors only).
    % Extend the structure for storing primary and secondary MBFs as follows:
    % Structure to be followed:   (MBF{1:Ndom},PrimIndx,DomainInd,solNum)
    % See FEKDDM-3.7: The following has to be changed to account for
    % multiple primary MBFs (as would be the case for multiple ports per
    % base domain) - then change the second dimension.
    % 2018.06.04: Increase now the size of PrimIsol to the global MoM size (not per domain any more).
    % 2018.06.07: The number of primaries / domain will not be 1 when we have interconnected domains
    % with generating sub-arrays (radiating case).
    max_primaries_per_domain = 1;
    if (~Solver_setup.disconnected_domains)
        % Number of primaries / domain - in the case of generating sub-arrays will depend on the location
        % of the element in the array. As a maximum, just take the maximum number of domains in any 
        % generating sub-array configuration. Loop over all sub-arrays and extract the max domain size
        max_subarray_domains = 0;
        for ii = 1:numGeneratingSubarrays
            max_subarray_domains = max(max_subarray_domains, ...
                length(Solver_setup.generating_subarrays.domains{ii}));            
        end%for ii=1:Solver_setup.generating_subarrays.number_of_domains
        max_primaries_per_domain = max_subarray_domains;

        % We should also allocate some space to store the Primaries associated with 
        % the generating sub-arrays. Note: we do not bother with a per Solution configuration
        % here, as the manner in which we generate these MBFs should make them useful for any
        % excitation configuration.
        mbfs.GenSubArrPrimIsol = complex(zeros(Nmom,max_primaries_per_domain,numGeneratingSubarrays,1));
        % Number of Generating Prim. MBFs / subarray
        mbfs.numGenSubArrPrimMBFs = zeros(numGeneratingSubarrays,1); 
    end%if

    mbfs.PrimIsol = complex(zeros(Nmom,max_primaries_per_domain,numArrayEls,numSols));    

    mbfs.numPrimMBFs = zeros(numArrayEls,numSols); % Number of Prim. MBFs / solution config.
    mbfs.numSecMBFs = zeros(numArrayEls,numSols);  % Number of Sec.  MBFs / solution config.
    if (Const.calcSecMBFs)
        % TO-DO: Adjust the second dimension of the secondary CBFs according
        % to the number that will be included - here, all neighbouring elements 
        % are accounted for (i.e. numArrayEls - 1 secondary MBFs will be
        % induced) - see also FEKDDDM-2 for additional features that are planned
        %                           (Ndom,  SecIndx    ,  DomainInd)
        % See issue FEKDDM-3.5 (and FEKDDM-10): Added now support for
        % multiple solution configurations. Threfore extend the structure
        % as follows:               (Ndom,  SecIndx    ,  DomainInd, # Sol. Configurations)
        % 2018.06.04: Increase now the size of SecIsol to the global MoM size (not per domain any more).
        % 2018.06.07: For connected arrays, we follow the same reasoning as above - i.e. the number of 
        % secondaries depend on the number of primaries that again depend on the number of generating sub-arrays.
        mbfs.SecIsol = complex(zeros(Nmom,max_primaries_per_domain*(numArrayEls-1),numArrayEls,numSols));
    end%if
    

    % 2018.06.03: We also support now interconnected domains. Only therefore 
    % pre-allocate certain MBF datastructures if we have a disjoint array problem.
    % Otherwise, these have to be calculated per domain.
    if (Solver_setup.disconnected_domains)

        % We are working with identical domains. Generate the LU decomposition
        % of the static interaction matrix of domain 1 beforehand and reuse
        % this in the following calculations
        % TO-DO: Danie, save this to the Temp directory and load it here to save time
        % See issue FEKDDM-6.2: Improved now the basis function numbering for
        % the array domain (work with bottom and top basis function offsets)
        % TO-DO: Assumed here are 1-to-1 mapping, i.e. Const.arrayMappingVector
        %        not yet used.        
        domain_indices = Solver_setup.rwg_basis_functions_domains{1};    
        % TO-DO: Actually use the calcZmn function here with the correct frequency index.
        [L,U] = lu(zMatrices.values(domain_indices, domain_indices));
    else

        message_fc(Const,sprintf('Generating sub-array primary MBFs'));

        % We have a generating sub-array here. We need to generate now the correct primaries
        % for each sub-array following the method as explained in [2].
        subarray_window = complex(zeros(Nmom, 1));        
        
        %  Loop over the number of sub-arrays
        for ii = 1:numGeneratingSubarrays
            
            % Extract the number of domains within this sub-array
            subarray_domains = Solver_setup.generating_subarrays.domains{ii};
            num_domains = length(subarray_domains);

            % Extract the unknowns associated with the entire sub-array (will be used below).
            subarray_basis_functions = ...
                Solver_setup.generating_subarrays.rwg_basis_functions_domains{ii};

            % Calculate the windowing function for this sub-array. Distinguish between 3 types:
            %   Sub-array 1 : Left edge sub-array
            %   Sub-array 2 : Centre sub-array
            %   Sub-array 3 : Right edge sub-array
            
            % Reset the window first.
            subarray_window(:) = 0.0;

            switch ii
                case 1
                    fprintf("  Sub-array (Type 1) : Left edge sub-array\n");
                    internal_domain = subarray_domains(1);
                case 2
                    fprintf("  Sub-array (Type 2) : Centre sub-array\n");
                    internal_domain = subarray_domains(2);
                case 3
                    fprintf("  Sub-array (Type 3) : Right edge sub-array\n");
                    internal_domain = subarray_domains(2);
                otherwise
                    % Error : we can only process 3 types of generating sub-arrays thus far
                    message_fc(Const,"Subarray type not supported");
                    error(['Subarray type not supported']);
            end

            % Extract the internal basis functions for the subarray, i.e. the one to which
            % the primary MBF will be mapped.
            domain_basis_functions_internal = ...
                Solver_setup.rwg_basis_functions_internal_domains{internal_domain};

            % Determine the basis functions on the overlapping region (i.e. on the interface 
            % of this element the surrounding elements)
            domain_basis_functions_interface = setdiff(...
                Solver_setup.rwg_basis_functions_domains{internal_domain}, ...
                Solver_setup.rwg_basis_functions_internal_domains{internal_domain});

            % Determine the basis functions external to this element (and also not on the interface)
            domain_basis_functions_external = setdiff(...
                subarray_basis_functions, ...
                Solver_setup.rwg_basis_functions_domains{internal_domain});

            subarray_window(domain_basis_functions_internal)  = 1.0;
            subarray_window(domain_basis_functions_interface) = 0.5;
            subarray_window(domain_basis_functions_external)  = 0.0;
            
            % We will be needing a vector to store the excited array
            % values
            yVectors_genPrim = complex(zeros(Nmom,1));

            % For each of the domains, we calculate primary MBFs by looping over the elements
            % and exciting each individually.
            count = 0;
            for m=1:num_domains % within the sub-array
                % Extract domain index from sub-array list
                domain_index = subarray_domains(m);

                % Extract the basis functions associated with the element internal to the 
                % domain (within the generating sub-array). This is needed to ensure that we excite only
                % this element, and also below to apply the correct windowing. 
                domain_basis_functions_internal = Solver_setup.rwg_basis_functions_internal_domains{domain_index};

                % Now, if we have a generating sub-array, i.e. for connected domains, then we excite ONLY the current element.                
                yVectors_genPrim(domain_basis_functions_internal) = yVectors.values(domain_basis_functions_internal,1);

                % Now let's calculate the primary MBF that is a result of exciting this element within the 
                % sub-array. TO-DO: Replace later with calcZmn.m
                [L,U] = lu(zMatrices.values(subarray_basis_functions, subarray_basis_functions, 1));

                b = L\yVectors_genPrim(subarray_basis_functions);
                
                count = count + 1;
                mbfs.GenSubArrPrimIsol(subarray_basis_functions,count,ii,1) = U\b; % U, already calculated above

                % Window now this MBF by applying the windowing function that was calculated above
                mbfs.GenSubArrPrimIsol(:,count,ii,1) = subarray_window.*mbfs.GenSubArrPrimIsol(:,count,ii,1);

            end%for m=1:num_domains        
            % Store the number of primary MBFs that were generated for this sub-array    
            mbfs.numGenSubArrPrimMBFs(ii) = count;

        end % ii = 1:numGeneratingSubarrays
    end

    % See issue FEKDDM-10: We added now support for multiple solution configurations.
    % Each solution configuration get its own set of primary and secondary MBFs, 
    % depending on the solution excitation configuration.
    for solNum = 1:numSols
    
        % Start timing (per solution)
        tic
        
        % ======================================
        % Setup the primary MBFs
        % ======================================
        
        % Generate the primary MBF: Jprim = (Zself)^(-1) * Vself
        mbfs.numPrimMBFs(:,solNum) = 0;
        
        for m=1:numArrayEls
            % We only generate a primary MBF if the array element is active
            if (Const.is_array_element_active(m,solNum))

                % Back-wards substitution with the part of the excitation vector local to this domain.
                % Note: if we have an interconnected domain problem, then this represents
                % the extended domain's solution.
                domain_basis_functions = Solver_setup.rwg_basis_functions_domains{m};

                % 2018.06.03: We also support now interconnected domains. Calculate
                % the LU decomposition here for the particular domain (will not be
                % the same for each element)
                if (~Solver_setup.disconnected_domains)
                    % TO-DO: Actually use the calcZmn function here with the correct frequency index.
                    [L,U] = lu(zMatrices.values(domain_basis_functions, domain_basis_functions));
                end

                mbfs.numPrimMBFs(m,solNum) = mbfs.numPrimMBFs(m,solNum) + 1;

                b = L\yVectors.values(domain_basis_functions,solNum);
                mbfs.PrimIsol(domain_basis_functions,1,m,solNum) = U\b; % U, already calculated above

                % Before we continue, we need to window the primary MBF here, if we are working with interconnected
                % domains:
                if (~Solver_setup.disconnected_domains && true)
                    % -- Windowing  (only if we have interconnected domains)
                    
                    % Let's first determine the BFs on the interface esssentially the difference between the 
                    % unknowns internal to the domain and that on the interface.
                    interface_basis_functions = setdiff(domain_basis_functions, ...
                        Solver_setup.rwg_basis_functions_internal_domains{m});

                    % Apply now windowing : Factor of a half.
                    mbfs.PrimIsol(interface_basis_functions,1,m,solNum) = 0.5.*mbfs.PrimIsol(interface_basis_functions,1,m,solNum);
                end%if
            end%if
        end%for
        
        % End timing
        mbfs.primGenTime(solNum) = toc;
        
        % ======================================
        % Setup the secondary MBFs (if included)
        % ======================================
        if (Const.calcSecMBFs)
            
            tic % Start timing
            
            % Generate the secondary MBF: Jsec = (Zself)^(-1) * Zcoupl * Jprim
            % for domain m, that is excited by the primary MBF on domain n
            for m=1:numArrayEls

                % Extract basis function indices of domain m.
                domain_m_basis_functions = Solver_setup.rwg_basis_functions_domains{m};

                % 2018.06.03: We also support now interconnected domains. Calculate
                % the LU decomposition here for the particular domain (will not be
                % the same for each element)
                if (~Solver_setup.disconnected_domains)
                    % TO-DO: Actually use the calcZmn function here with the correct frequency index.
                    [L,U] = lu(zMatrices.values(domain_m_basis_functions, domain_m_basis_functions));

                    % In addition, we will also need to apply a windowing to the secondary MBF - similar to
                    % what was done for the primary MBF. This will be done below. Calculate here first the
                    % rwg coefficients on the interface (overlapping region) so that we can apply a windowing
                    % function.
                    if (true)
                        % -- Windowing  (only if we have interconnected domains)                        
                        % Let's first determine the BFs on the interface esssentially the difference between the 
                        % unknowns internal to the domain and that on the interface.
                        interface_basis_functions = setdiff(domain_m_basis_functions, ...
                            Solver_setup.rwg_basis_functions_internal_domains{m});                        
                    end%if
                end

                % Back-wards substitution with the part of the excitation vector
                % local to this domain
                count = 0;
                for n=1:numArrayEls 
                    if (m == n)
                        %ignore self-coupling - primary MBF already calculated
                        continue;
                    end%if

                    % Keep track of the number of secondary MBFs on domain
                    % m, if domain n has a primary to excite such a secondary MBF
                    if (Const.is_array_element_active(n,solNum))
                        count = count + 1;
                        mbfs.numSecMBFs(m,solNum) = count;                                  

                        % Extract basis function indices of domain n. See [1], if we have interconnected domains, 
                        % then we use the internal RWGs (i.e. not including the unknowns on the interface) for
                        % the source/basis MBFs.
                        
                        if (~Solver_setup.disconnected_domains)
                            domain_n_basis_functions = Solver_setup.rwg_basis_functions_internal_domains{n};
                        else
                            domain_n_basis_functions = Solver_setup.rwg_basis_functions_domains{n};
                        end                        

                        % Calculate the coupling matrix
                        % TO-DO: Actually use the calcZmn function here with the correct frequency index.
                        Zcoupl = zMatrices.values(domain_m_basis_functions,domain_n_basis_functions);
                        % Calculate the field coupling to domain m, using primary MBF from domain n
                        % Note: cf Eq. (4) in [1], the excitation vector resulting from the mutual 
                        % coupling has a negative sign! Be careful here when reusing these secondary 
                        % MBFs in the Jacobi Iterative Solver - then no negative sign is used.                        
                        Vcoupl =  - Zcoupl * mbfs.PrimIsol(domain_n_basis_functions,1,n,solNum);

                        % Solve now for the secondary induced MBF using the previously calculated 
                        % LU-decomposition of Zmm (stored in L and U)
                        b = L\Vcoupl;
                        mbfs.SecIsol(domain_m_basis_functions,count,m,solNum) = U\b;

                        % Window the secondary MBF here, if we are working with interconnected domains:
                        if (~Solver_setup.disconnected_domains && true)
                            % Apply now windowing : Factor of a half.
                            mbfs.SecIsol(interface_basis_functions,count,m,solNum) = 0.5.*mbfs.SecIsol(interface_basis_functions,count,m,solNum);
                        end%if
                    end%if
                end%for
            end        
            mbfs.secGenTime(solNum) = toc; % End timing
        else
            mbfs.secGenTime(solNum) = 0.0;
        end%if
    
% --------------------------------------------------------------------------------------------------
        % 2015-08-18: Reduce and orthonormalize MBFs (taken from CEASER2p5)
        if (Const.useMBFreduction)
            tic
            message_fc(Const,sprintf('Reduce and orthonormalise MBFs'));
            for m=1:numArrayEls            
                % Put all the MBFs in a column augmented matrix
                if (Const.calcSecMBFs)
                    origMBFs = [mbfs.PrimIsol(:,1,m,solNum) mbfs.SecIsol(:,1:mbfs.numSecMBFs(m,solNum),m,solNum)];
                else
                    origMBFs = mbfs.PrimIsol(:,1,m,solNum);
                end
                if (Const.debug)
                    message_fc(Const,['Number of initially generated CBFs: ' num2str(size(origMBFs,2))]);
                    message_fc(Const,'Reducing this number based on the user specified threshold...');
                end%if

                fcdString = sprintf('cbfm');
                redMBFs = reduceMBFset(Const, origMBFs, fcdString);
                %redMBFs = reduceMBFset(origMBFs,Const.MBFthreshold,Const.MBFplotSVspectrum);
                mbfs.RedIsol(:,1:size(redMBFs,2),m,solNum) = redMBFs;
                mbfs.numRedMBFs(m,solNum) = size(redMBFs,2);
                if (Const.debug)
                    message_fc(Const,['Number of retained orthonormal CBFs: ' num2str(size(redMBFs,2))]);
                end%if                
            end%for
            mbfs.svdTime(solNum) = toc; % End timing
        else
            mbfs.svdTime(solNum) = 0.0;
            message_fc(Const,sprintf('  No MBF reduction (+ orthonormalisation)'));
        end %if (Const.useMBFreduction)
% --------------------------------------------------------------------------------------------------
    end %for solNum = 1:numSols

    message_fc(Const,sprintf('Finished MBF generator'));
    % Output per solution configuration data
    totprimGenTime = 0;
    totsecGenTime  = 0;
    totsvdTime = 0;
    for solNum = 1:numSols
        % Calculate the total number of Primary and Secondary MBFs for the solution
        totPrimMBFs = 0;
        totSecMBFs  = 0;
        mbfs.totRedMBFs  = 0; % Store this for later reuse
        for n=1:numArrayEls
            totPrimMBFs = totPrimMBFs + mbfs.numPrimMBFs(n,solNum);
            totSecMBFs  = totSecMBFs  + mbfs.numSecMBFs(n,solNum);
            if (Const.useMBFreduction)
                mbfs.totRedMBFs = mbfs.totRedMBFs + mbfs.numRedMBFs(n,solNum);
            end%if
        end
        message_fc(Const,sprintf('Total number of Primary MBFs %d , Number of Second. MBFs %d for Sol. %d of %d',...
            totPrimMBFs, totSecMBFs, solNum, numSols));
        if (Const.useMBFreduction)
            message_fc(Const,sprintf('Total number of Reduced MBFs %d for Sol. %d of %d',...
            mbfs.totRedMBFs, solNum, numSols));
        end%if
        message_fc(Const,sprintf('Times for calculating Primary MBFs %f sec. and Second. MBFs %f sec. for Sol. %d of %d',...
            mbfs.primGenTime(solNum), mbfs.secGenTime(solNum), solNum, numSols));
        if (Const.useMBFreduction)
            message_fc(Const,sprintf('Times for calculating reduced MBFs (SVD) %f sec. for Sol. %d of %d',...
                mbfs.svdTime(solNum), solNum, numSols));
        end%if
        % Calculate the total time
        totprimGenTime = totprimGenTime + mbfs.primGenTime(solNum);
        totsecGenTime  = totsecGenTime  + mbfs.secGenTime(solNum);
        totsvdTime  = totsvdTime  + mbfs.svdTime(solNum);
    end
    % Output total time
    message_fc(Const,sprintf('Total Times for calculating Primary MBFs %f sec. and Second. MBFs %f sec.',totprimGenTime, totsecGenTime));
    if (Const.useMBFreduction)
        message_fc(Const,sprintf('Total Times for reducing MBFs (SVD) %f sec. ',totsvdTime));
    end%if    

    % We need to store also the MBF total time here - reported at the end of the CBFM solver when results are written
    % to file
    mbfs.totTime = totprimGenTime + totsecGenTime + totsvdTime;