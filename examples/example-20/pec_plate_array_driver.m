% Author: Danie Ludick (dludick@sun.ac.za)
% Project: PEC plate array example
%
% Note: Each project directory / example directory needs to have a sunem_initialise.m
% script, that is used to setup the correct environment for that example.
%
% Refer to the /doc folder for more information

% --------------------------------------------------------------------------------------------------
% Initialise the environment
% --------------------------------------------------------------------------------------------------
% Project output directory: './dipoles/'
% Debug: True/False
Const = sunem_initialise('pec_plate_array',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver       = true;
Const.runCBFMsolver      = false;
Const.runJacobisolver    = false;
Const.runIFBMoMsolver    = true;
Const.runDGFMsolver      = false;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'pec_plate_array.mat';
Const.FEKOstrfilename          = 'pec_plate_array.str';
Const.FEKOrhsfilename          = 'pec_plate_array.rhs';
Const.FEKOoutfilename          = 'pec_plate_array.out';
Const.FEKOefefilename          = 'pec_plate_array.efe';
Const.FEKOffefilename          = 'pec_plate_array.ffe';

% --------------------------------------------------------------------------------------------------
% Define output files for transferring expansion coefficients back to FEKO data
% --------------------------------------------------------------------------------------------------
Const.SUNEMifbmomstrfilename   = '';%'ifbmom_pec_plate_array.str';
Const.SUNEMcbfmstrfilename     = '';

% --------------------------------------------------------------------------------------------------
% Define additional program flow constants
% --------------------------------------------------------------------------------------------------
% TO-DO: Setup some documentation for this - also assign default values if they are not
% defined explicitely here.
Const.no_mutual_coupling_array = false; % Deactivate coupling between domains.
Const.calcSecMBFs = false;     % For MBF based solvers
Const.useMBFreduction = true;  % SVD applied after the MBFs are generated to retain an orthonormal set
Const.MBFthreshold = -1;       % Threshold used for the SVD reduction of the MBFs
Const.IFBalg = 7;              % Jacobi iterations (7). Adaptive MBF (14).
Const.IFB_iterations = 3;     % Number of Jacobi iterations. (TO-DO: Ellaborate special meaning, e.g. -1)
                               % which then looks at Const.IFB_convergence_threshold_percentage;
Const.IFB_convergence_threshold_percentage = 1E-3;                                
Const.IFB_CBFs = -1;           % TO-DO: Recheck this - essentially for the Adaptive MBF the number of MBFs
                               % to use during each iteration
Const.IFB_debug = 1;
Const.cache_Z0_V0 = false;     % Precompute the Z0 and V0 terms
Const.use_DGFM_start = false;  % Use the DGFM to calculate the initial (0th) solution

% --------------------------------------------------------------------------------------------------
% Read the MoM matrix equation from the file
% --------------------------------------------------------------------------------------------------
[Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const);

% --------------------------------------------------------------------------------------------------
% Parse the setup files to extract the frequency sweep, the geometry and basis function setup 
% --------------------------------------------------------------------------------------------------
% TO-DO: At a later stage we can also add other meshing / geometry
% preprocessxing, e.g. Gmsh or GiD. For now the solver setup is read from FEKO.
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectors);

% --------------------------------------------------------------------------------------------------
% Run the EM solver
% --------------------------------------------------------------------------------------------------
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);