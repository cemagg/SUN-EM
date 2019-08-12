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

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'strip_dipoles_random_array.mat';
Const.FEKOstrfilename          = 'strip_dipoles_random_array.str';
Const.FEKOrhsfilename          = 'strip_dipoles_random_array.rhs';
Const.FEKOoutfilename          = 'strip_dipoles_random_array.out';
Const.FEKOefefilename          = 'strip_dipoles_random_array.efe';
Const.FEKOffefilename          = 'strip_dipoles_random_array.ffe';

% --------------------------------------------------------------------------------------------------
% Define output files for transferring expansion coefficients back to FEKO data
% --------------------------------------------------------------------------------------------------
Const.SUNEMmomstrfilename      =  'SUNEM_MOM_strip_dipoles_random_array.str';%'sunem_mom_bow_tie_array.str';

% --------------------------------------------------------------------------------------------------
% Define additional program flow constants
% --------------------------------------------------------------------------------------------------
Const.no_mutual_coupling_array = false; % Deactivate coupling between domains.
Const.calcSecMBFs = true;               % For MBF based solvers
Const.useMBFreduction = true;           % SVD applied after the MBFs are generated to retain an orthonormal set
Const.MBFthreshold = 1000;              % Threshold used for the SVD reduction of the MBFs

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
% (Note: We either pass our own (internal) matrices, or those read from FEKO). For this particular
% array configuration, we are not yet supporting radiating elements. But as we are consistent with the
% FEKO mesh, we can just use the FEKO RHS vector.
[Solution] = runEMsolvers(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO);

