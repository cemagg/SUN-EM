%function [Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors)
function [SolutionFEKO] = runEMsolvers(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO)
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
    %           MoM, CBFM
    %
    %   Description:
    %       Runs various EM solvers based on the data that was read / parsed.
    %       from the FEKO *.out, *.mat, *.str and *.rhs files (for now).

    narginchk(5,5);
    %narginchk(6,6);
  
    % Initialise the return values
 %Solution  = [];
 SolutionFEKO = [];

%     % -- MoM    
%     if (Const.runMoMsolver)        
%    Solution.mom = runMoMsolver(Const, Solver_setup, yVectors, zMatrices, xVectors);
%     end%if
    % -- MoM    
    if (Const.runMoMsolver)        
         SolutionFEKO.mom = runMoMsolver(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO);
    end%if



    % -- C++ MoM (TO-DO: Tameez, perhaps this is a better spot to call your entire MoM C++ solver)
     %if (Const.runMoMsolver)        
    %     Solution.mom = runCPPMoMsolver(Const, Solver_setup, zMatrices, yVectors, xVectors);
   %  end%if    
 
    % -- CBFM
    if (Const.runCBFMsolver)
        % First generate the MBFs (separate routine, as we might be using
        % different MBF generation techniques)
        [mbfs] = runMBFgenerator(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO);
        
        % Now reduced matrix setup + solution
        Solution.cbfm = runCBFMsolver(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO, mbfs);
    end%if

    % -- Jacobi solver (only two iterations, as used by the i-DGFM)
    if (Const.runJacobisolver)
        Solution.jack = runJACKITsolver(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO, mbfs);
    end%if

    % -- IFB MoM solver (similar to Jacobi)
    if (Const.runIFBMoMsolver)
        Solution.ifbmom = runIFBMoMsolver(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO);
    end%if
    
    % -- FPGA-based Iterative Jacobi solver (MEng of Caleb Mnisi - 2020)
    if (Const.runFPGAJacobisolver)
        Solution.fpgajacobi = runFPGAJacobisolver(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO);
    end%if
    
    % -- DGFM solver
    if (Const.runDGFMsolver)
        ngf = []; % TO-DO: Add back NGF once supported.
        % Depending on the DGFM weighting coefficients, we might need to
        % run the MBF generator
        if (Const.DGFMweightVectorCalcScheme == 3)
            [mbfs] = runMBFgenerator(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO);
        else
            mbfs = [];
        end

        SolutionFEKO.dgfm = runDGFMsolver(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO, mbfs, ngf);  
        
        if (Const.useDGFMinterpolation)
            % In order to extrapolate the DGFM obtained MBF coefficients, we
            % first need to visualise them on the array processing grid to
            % see whether they are smooth. Some initial tests can then also be
            % done here.
            SolutionFEKO.dgfm = dgfmMBFinterpolate(Const, Solver_setup, yVectorsFEKO, SolutionFEKO.dgfm, mbfs);
        end %if (Const.useDGFMinterpolation)
        
    end%if

    message_fc(Const,sprintf('Finished EM Solvers'));