function [Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors)
    %runEMsolvers
    %   Usage:
    %       [Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors)] = runEMsolvers(Const)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containging settings of which solver to run
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details
    %       zMatrices
    %           The Z-matrices data OR can be empty
    %       yVectors
    %           The Yrhs-vector data OR can be empty
    %       xVectors
    %           The Xsol-vector data (i.e. MoM solution of FEKO) OR can be empty
    %
    %   Output Arguments:
    %      Solution
    %           Structs containing solutions and timing data of EM solvers, e.g.
    %           MoM, HARP
    %
    %   Description:
    %       Runs various EM solvers based on the data that was read / parsed.
    %       from the FEKO *.out, *.mat, *.str and *.rhs files (for now).

    error(nargchk(5,5,nargin));

    % Initialise the return values
    Solution  = [];

    % -- MoM (FEKO specific)
    % This solver actually extracts the needed parameters from the zMatrices and yVectors structs.
    % TO-DO: We need to think about how to refactor to include also our local GMoM solver
    if (Const.runMoMsolver)        
        mom = runMoMsolver(Const, zMatrices, yVectors, xVectors);
        Solution.mom = mom;
    end%if

    % -- HARP
    if (Const.runHARPsolver)
        harp = runHARPsolver(Const, Solver_setup, zMatrices, yVectors, xVectors);
        Solution.harp = harp;
    end%if

    message_fc(Const,sprintf('Finished EM Solvers'));