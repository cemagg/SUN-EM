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
Const = sunem_initialise('dipole_antenna_array',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver          = true;
Const.runCBFMsolver         = true;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'dipole_antenna_array.mat';
Const.FEKOstrfilename          = 'dipole_antenna_array.str';
Const.FEKOrhsfilename          = 'dipole_antenna_array.rhs';
Const.FEKOoutfilename          = 'dipole_antenna_array.out';
%Const.FEKOefefilename          = 'dipole_antenna_array.mat.efe';
%Const.FEKOffefilename          = 'dipole_antenna_array.mat.ffe';

% Added also array layout specification
% Const.arrayLayoutfilename      = 'array_layout.xml';

% --------------------------------------------------------------------------------------------------
% Define output files for transferring expansion coefficients back to FEKO data
% --------------------------------------------------------------------------------------------------
Const.SUNEMmomstrfilename      =  '';%'dipole_antenna_array.str';

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

