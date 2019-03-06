% Author: Danie Ludick (dludick@sun.ac.za)
% Project: PEC plate example
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
Const = sunem_initialise('random_strip_dipoles',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver          = true;
Const.runCBFMsolver         = true;
Const.runJacobisolver       = false;
Const.runIFBMoMsolver       = false;
Const.runDGFMsolver         = true;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'strip_dipoles_random_array.mat';
Const.FEKOstrfilename          = 'strip_dipoles_random_array.str';%'strip_dipoles_random_array.str';
Const.FEKOrhsfilename          = 'strip_dipoles_random_array.rhs';
Const.FEKOoutfilename          = 'strip_dipoles_random_array.out';
Const.FEKOefefilename          = 'strip_dipoles_random_array.efe';
Const.FEKOffefilename          = 'strip_dipoles_random_array.ffe';

% Added also array layout specification
Const.arrayLayoutfilename      = 'array_layout.xml';

% --------------------------------------------------------------------------------------------------
% Define output files for transferring expansion coefficients back to FEKO data
% --------------------------------------------------------------------------------------------------
Const.SUNEMmomstrfilename      =  'SUNEM_MOM_strip_dipoles_random_array.str';%'sunem_mom_bow_tie_array.str';
Const.SUNEMdgfmstrfilename     =  'DGFM_strip_dipoles_random_array.str';%'sunem_dgfm_bow_tie_array.str';

% --------------------------------------------------------------------------------------------------
% Define additional program flow constants
% --------------------------------------------------------------------------------------------------
% TO-DO: Setup some documentation for this
Const.no_mutual_coupling_array = false; % Deactivate coupling between domains.
Const.calcSecMBFs = true;      % For MBF based solvers
Const.useMBFreduction = true;   % SVD applied after the MBFs are generated to retain an orthonormal set
Const.MBFthreshold = 1000;      % Threshold used for the SVD reduction of the MBFs (typicall 1000)
Const.IFBalg = 14;              % Jacobi iterations (7). Adaptive MBF (14).
Const.IFB_iterations = 10;      % Number of Jacobi iterations. (TO-DO: Ellaborate special meaning, e.g. -1)
                                % which then looks at Const.IFB_convergence_threshold_percentage;
Const.IFB_convergence_threshold_percentage = 1E-3;                                
Const.IFB_CBFs = -1;            % TO-DO: Recheck this - essentially for the Adaptive MBF the number of MBFs
                                % to use during each iteration
Const.IFB_debug = 1;
Const.cache_Z0_V0 = false;      % Precompute the Z0 and V0 terms
Const.use_DGFM_start = false;   % Use the DGFM to calculate the initial (0th) solution

Const.useEDM = true;            % Use the Equivalent Dipole Method (EDM) to accelerate the MoM Z-matrix
                                % calculation

% ------------------
%  DGFM coefficients
% ------------------
Const.useDGFMmethod = 1;                 % 1: Local matrices. 2: Global matrices (not supported anymore)
Const.storeZact     = 1;
Const.DGFMweightVectorCalcScheme = 3;    % 0: Uniform coefficients / array excitation law. 
                                         % 1: Ratio of the applied excitation coefficients (TO-DO: check)
                                         % 2: Use Isol reference. (Great for testing)
                                         % 3: Use the Jacobi-Generated CBFs (see [DL2013], Appendix A)
                                         
Const.useDGFMinterpolation = 1;          % 0: No interpolation
                                         % 1: Use interpolation (MATLAB internal scattered interpolant method)
                                         % 2: TO-DO: Krigin interpolation
                                         
Const.DGFMinterpolationSamplingFactor = 0.2; % The sampling factor of the array for the interpolation
                                             % where 1 corresponds to a fully sampled array.
                                             % Each sample corresponds to a full DGFM calculation.
                                         
% --------------------------------------------------------------------------------------------------
% Read the MoM matrix equation from the file
% --------------------------------------------------------------------------------------------------
[Const, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO] = extractFEKOMoMmatrixEq(Const);

% --------------------------------------------------------------------------------------------------
% Parse the setup files to extract the frequency sweep, the geometry and basis function setup 
% --------------------------------------------------------------------------------------------------
% TO-DO: At a later stage we can also add other meshing / geometry
% preprocessxing, e.g. Gmsh or GiD. For now the solver setup is read from FEKO.
[Const, Solver_setup] = parseFEKOoutfile(Const, yVectorsFEKO);

% --------------------------------------------------------------------------------------------------
% Run the EM solver 
% --------------------------------------------------------------------------------------------------
[Solution] = runEMsolvers(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO);


