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
        message_fc(Const,sprintf('(*** Mutual coupling between array elements ignored ***)'));
        % Simple - just switch off calculating secondary MBFs
        Const.calcSecMBFs = false;
    else
        % Calculate also secondary CBFS
        Const.calcSecMBFs = true;
    end%if
        
    % Initialise the return values
    mbfs  = [];    
    Nmom        = Solver_setup.num_mom_basis_functions;    % Total number of basis functions for whole problem
    Nngf        = Solver_setup.num_ngf_basis_functions;    % Number of basis functions for NGF domain
    Ndom        = Solver_setup.mom_basis_functions_per_array_element;   % Number of basis functions per array element
    numArrayEls = Solver_setup.num_finite_array_elements;  % The number of array elements
    
    numSols     = xVectors.numSols;             % The number of solutions configurations
    
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
    mbfs.PrimIsol = complex(zeros(Ndom,1,numArrayEls,numSols));
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
        mbfs.SecIsol = complex(zeros(Ndom,numArrayEls-1,numArrayEls,numSols));
    end%if
    
    % If we are reducing and orthonormalising the MBFs by using the SVD,
    % then store these in a differents array than mbfs.PrimIsol and
    % mbfs.SecIsol
    if (Const.useMBFreduction)
        % Note, we are unsure at this stage how many reduced MBFs we will
        % retain. Just allocated enough space here (numArrayEls should do
        % the trick).
        mbfs.RedIsol = complex(zeros(Ndom,numArrayEls,numArrayEls,numSols));
        mbfs.numRedMBFs = zeros(numArrayEls,numSols); % Number of Red. MBFs / array element / solution config.
    end%if
        
    % We are working with identical domains. Generate the LU decomposition
    % of the static interaction matrix of domain 1 beforehand and reuse
    % this in the following calculations
    % TO-DO: Danie, save this to the Temp directory and load it here to save time
    % See issue FEKDDM-6.2: Improved now the basis function numbering for
    % the array domain (work with bottom and top basis function offsets)
    % TO-DO: Assumed here are 1-to-1 mapping, i.e. Const.arrayMappingVector
    %        not yet used.
    %domain_bot = Const.arrayElBasisFunctRange(1,1);
    %domain_top = Const.arrayElBasisFunctRange(1,2);
    %[L,U] = lu(zMatrices.values(domain_bot:domain_top,domain_bot:domain_top));
    
    domain_indices = Solver_setup.rwg_basis_functions_domains{1};
    
    [L,U] = lu(zMatrices.values(domain_indices, domain_indices) );
    
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
                mbfs.numPrimMBFs(m,solNum) = mbfs.numPrimMBFs(m,solNum) + 1;
                
                % Back-wards substitution with the part of the excitation vector local to this domain                
                domain_indices = Solver_setup.rwg_basis_functions_domains{m};
                
                b = L\yVectors.values(domain_indices,solNum);
                mbfs.PrimIsol(:,1,m,solNum) = U\b; % U, already calculated above
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
                    if (Const.isArrayElementActive(n,solNum))
                        count = count + 1;
                        mbfs.numSecMBFs(m,solNum) = count;
                        % Calculate the coupling matrix
                        domain_bot_m = Const.arrayElBasisFunctRange(m,1);
                        domain_top_m = Const.arrayElBasisFunctRange(m,2);
                        domain_bot_n = Const.arrayElBasisFunctRange(n,1);                
                        domain_top_n = Const.arrayElBasisFunctRange(n,2);
                        Zcoupl = zMatrices.values(domain_bot_m:domain_top_m,domain_bot_n:domain_top_n);
                        % Calculate the field coupling to domain m, using primary MBF from domain n
                        % Note: cf Eq. (4) in [1], the excitation vector resulting from the mutual 
                        % coupling has a negative sign! Be careful here when reusing these secondary 
                        % MBFs in the Jacobi Iterative Solver - then no negative sign is used.
                        Vcoupl =  - Zcoupl * mbfs.PrimIsol(:,1,n,solNum);
                        % Solve now for the secondary induced MBF using the previously calculated 
                        % LU-decomposition of Zmm (stored in L and U)
                        b = L\Vcoupl;
                        mbfs.SecIsol(:,count,m,solNum) = U\b;
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
                origMBFs = [mbfs.PrimIsol(:,1,m,solNum) mbfs.SecIsol(:,1:mbfs.numSecMBFs(m,solNum),m,solNum)];
                if (Const.debug)
                    message_fc(Const,['Number of initially generated CBFs: ' num2str(size(origMBFs,2))])
                    message_fc(Const,'Reducing this number based on the user specified threshold...')
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

