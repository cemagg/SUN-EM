function [jack] = runJACKITsolver(Const, Solver_setup, zMatrices, yVectors, xVectors, mbfs)
    %runJACKITsolver
    %   Date: 30.11.2013
    %   Usage:
    %       [jack] = runJACKITsolver(Const, Solver_setup, zMatrices, yVectors, xVectors)
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
    %       jack
    %           Struct containing the Jacobi Iterations solution and timing data
    %
    %   Description:
    %       Calculates the Jacobi Iterative Solver, as described in [1] and
    %       [2] and adapted for the DGFM in [3].
    %
    %   Assumptions:
    %        - All domains are the same size (i.e. contains the same number of unknowns)
    %
    %   TO-DO:
    %        - Currently only a single iteration is performed - we need to
    %          extend this to more than 1 iteration.
    %
    %   References:
    %   [1] Y. Brand, A. K. Skrivervik, J. R. Mosig, F. E. Gardoil, "New
    %       iterative integral equation technique for multilayered printed
    %       array antennas," in Mathematical Methods in Electromagnetic Theory,
    %       Kharkov, Ukraine, Jun, 1998, pp. 615-617.
    %   [2] Brand, Yan, Anja K. Skrivervik, and Juan R. Mosig. "An iterative scheme solution
    %       for the analysis of printed arrays." Microwave and Optical Technology Letters
    %       16.2 (1997): 106-115.
    %   [3] A. C. Polycarpou, "Evaluation of stationary block iterative
    %       techniques for the solution of finite arrays using fe-bi method
    %       and domain decomposition," in Proc. European Conference on
    %       Antennas and Propagation (EuCAP), Nice, France, Nov. 2006,
    %       pp. 1-6
    %   [4] D.J. Ludick, R. Maaskant, et. al, "Efficient Analysis of Large
    %       Irregular Antenna Arrays using the Domain Green's Function Method",
    %       Special Issue on Antennas Propagation Submission, 2013 (pending)
    %   =======================
    %   Written by Danie Ludick on June 24, 2013.
    %   Last updated on 2018.06.05.
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    narginchk(6,6);

    % Formating based on where the routine is called from
    if (Const.runJACKITfromDGFM)
        message_fc(Const,sprintf('  Running Jacobi Iterative solver'));
        if (Const.no_mutual_coupling_array)
            message_fc(Const,sprintf('  (*** Mutual coupling between array elements ignored ***)'));
        end%if
    else
        message_fc(Const,' ');
        message_fc(Const,'------------------------------------------------------------------------------------');
        message_fc(Const,sprintf('Running Jacobi Iterative solver'));
        if (Const.no_mutual_coupling_array)
            message_fc(Const,sprintf('(*** Mutual coupling between array elements ignored ***)'));
        end%if
    end%if

    if (Const.no_mutual_coupling_array)
        % If now mutual coupling should be accounted in the array environment, then just switch off calculating secondary MBFs
        Const.calcSecMBFs = false;
    end%if

    % Start timing
    tic

    % Initialise the return values
    jack  = [];
    jack.name = 'jack';
    Nmom        = Solver_setup.num_mom_basis_functions;  % Total number of basis functions for whole problem
    Nngf        = Solver_setup.num_ngf_basis_functions;  % Number of basis functions for NGF domain
    %Ndom        = Solver_setup.max_mom_basis_functions_per_array_element;  % Number of basis functions per array element
    numArrayEls = Solver_setup.num_finite_array_elements;   % The number of array elements

    numSols     = xVectors.numSols;             % The number of reference solutions
    jack.numSols = numSols;                     % Calculate a solution for each configuration
    jack.Isol   = complex(zeros(Nmom,numSols));   % The Jacobi method is used to calculate the array solution
    
    % Extract the solutions from which to start and end our process
    solStart = Const.solStart;
    solEnd   = Const.solEnd;

    % Note: We actually only perform a single Jacobi iteration here. Starting with the MBFs we generated
    % in runMBFgenerator. TO-DO: Add more iterations, or alternatively the CBFM for adaptive generation.
    jack.numIter = 1;
    jack.relError = zeros(numSols);
    
    % TO-DO: See FEKDDM-4.3: Check convergence of the Jacobi Iterative Solver
    if (checkJACKITconvergence(Const, zMatrices) ~= true)
       error (['No convergence possible for Jacobi Iter. Solver']);
       message_fc(Const,sprintf('No convergence possible for Jacobi Iter. Solver'));
    end

    % See issue FEKDDM-10: We added now support for multiple solution
    % configurations. Each of the solution configurations has given rise to
    % a suitable weightVectors (based on the excitation-law, exact current,
    % etc. for that specific solution). We need to repeat the DGFM active
    % impedance matrix calculation of each array element for each solution
    % configuration.
    for solNum = solStart:solEnd
        % -- Use the primary and secondary MBFs generated in the generic MBF routine
        for p=1:numArrayEls

            % Extract the basis functions for the domain
            % Note: if we have an interconnected domain problem, then this represents
            % the extended domain's solution.
            domain_basis_functions = Solver_setup.rwg_basis_functions_domains{p};

            % Add the contribution of the primary MBFs
            for prim_ii=1:mbfs.numPrimMBFs(p,solNum)
                % Note, assumed here for mbfs.PrimIsol, is that each of the domains only have 1
                % primary MBF. The last term in the following assignment has to be changed to:
                % mbfs.PrimIsol(:,m,p)
                jack.Isol(domain_basis_functions,solNum) = jack.Isol(domain_basis_functions,solNum) + ...
                    mbfs.PrimIsol(domain_basis_functions,prim_ii,p,solNum);
            end%for n=1:numArrayEls
            if (~Const.calcSecMBFs)
                % If there are no Secondary MBFs, then continue with the
                % next array element (only primary MBFs then taken into
                % account)
                continue;
            end%if

            %  Add the contribution of the secondary MBFs
            secIndx = 0;
            for sec_ii=1:mbfs.numSecMBFs(p,solNum)
                % See FEKDDM-4.1 and also comment in runMBFgenerator: The negative sign has already been
                % accounted for when the secondary MBFs are generated.
                jack.Isol(domain_basis_functions,solNum) = jack.Isol(domain_basis_functions,solNum) + mbfs.SecIsol(domain_basis_functions,sec_ii,p,solNum);
            end%for
        end%for

        % See issue FEKDDM-6.2: If the NGF-enhanced DGFM is used, then we
        % compare only the current that is on the dynamic domain, i.e. the
        % finite array. This routine might also be called from the IFB-DGFM solver (therefore used now Const.domAoffset instead of Nngf)
        jack.relError(solNum) = calculateErrorNormPercentage(xVectors.Isol(:,solNum), jack.Isol(:,solNum));
    end%for solNum = 1:numSols

    % End timing
    jack.solTime = toc;

    % Write the Jacobi iteration solution to a ASCII str file, so that it can be read
    % again by FEKO (for plotting in POSTFEKO) - only if requested (i.e. if the filename is defined)
    if (~isempty(Const.SUNEMjackstrfilename))
        writeSolToFile(Const, jack);
    end%if

    % Do output for each solution configuration
    for solNum = solStart:solEnd

        % Formating based on where the routine is called from
        if (Const.runJACKITfromDGFM)
            % Note, use here jack.numIter+1 because we include also the 0th iteration
            message_fc(Const,sprintf('  Finished Jacobi Iterative solver in %f sec. with %d iterations for Sol. %d of %d', ...
                jack.solTime, jack.numIter+1,solNum,numSols));
            message_fc(Const,sprintf('  Rel. error norm. compared to FEKO MoM sol. %f percent for Sol. %d of %d',...
                jack.relError(solNum),solNum,numSols));
        else
            message_fc(Const,sprintf('Finished Jacobi Iterative solver in %f sec. with %d iterations for Sol. %d of %d',...
                jack.solTime, jack.numIter+1,solNum,numSols));
            message_fc(Const,sprintf('Rel. error norm. compared to FEKO MoM sol. %f percent for Sol. %d of %d',...
                jack.relError(solNum),solNum,numSols));
        end%if
    end