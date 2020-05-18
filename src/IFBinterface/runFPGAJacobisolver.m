function [fpgajacobi] = runFPGAJacobisolver(Const, Solver_setup, zMatrices, yVectors, xVectors)
    %runFPGAJacobisolver
    %   Usage:
    %       [fpgajacobi] = runFPGAJacobisolver(Const, Solver_setup, zMatrices, yVectors, xVectors)
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
    %
    %   Output Arguments:
    %       fpgajacobi
    %           Structs containing the solution and timing data
    %
    %   Description:
    %       Runs the iterative Jacobi solver, with the option of having the FPGA enabled.
    %
    %   Notes:
    %       
    %   References: 
    %      [1] R. Maaskant, Ludick, D.J, Botha, M.M., D. Davidson, "Analysis of Finite Antennas Arrays
    %          using the CBFM enhanced Jacobi Method,", AWPL submission 2016 (work in progress)

    nargchk(5,5,nargin);

    % Some initialisations
    Nmom = Solver_setup.num_mom_basis_functions;                   % Total number of basis functions for whole problem
    numArrayElements = Solver_setup.num_finite_array_elements;     % The number of array elements

    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running FPGA Jacobi solver'));
    message_fc(Const,sprintf('Total basis functions   : %d.',Nmom));
    message_fc(Const,sprintf('Number of array elements: %d.',numArrayElements));
    message_fc(Const,sprintf('Number of Jacobi Iterations: %d.',Const.IFB_iterations));
    message_fc(Const,sprintf('Use FPGA : %d.',Const.useFPGA));

    % Set here a local debug flag on or off
    local_debug_flag = false;

    % Initialisations
    fpgajacobi  = [];
    fpgajacobi.name = 'fpgajacobi';

    Ntot  = Nmom;

    fpgajacobi.Isol = complex(zeros(Ntot,1));

    % Set the convergence threshold here (when Const.IFB_iterations = -1)
    eps_percent = Const.IFB_convergence_threshold_percentage;
    k_iter = Const.IFB_iterations;

    % Start timing
    fpgajacobi.solTime = 0;
    tic

    % ----------------------------------
    % TO-DO: Caleb, complete this routine .... (have a look at runfpgajacobisolver with Const.IFBalg equal to 7.
    % ----------------------------------

     % Stop pre-computation timing
    fpgajacobi.solTime = fpgajacobi.solTime + toc;


    % TO-DO: Can add more output here as needed.
    message_fc(Const,sprintf('Finished FPGAJacobisolver solver in %f sec.',fpgajacobi.totsolTime));