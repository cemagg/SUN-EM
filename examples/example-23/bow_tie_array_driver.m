% Author: Danie Ludick (dludick@sun.ac.za)
% Project: Bow Tie array simulation using MoM and DGFM
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
Const = sunem_initialise('bow_tie_antenna',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver       = true;
Const.runDGFMsolver      = false;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'bowtie-array_MoM.mat'; % Z-matrix calculated by FEKO
Const.FEKOstrfilename          = 'bowtie-array_MoM.str'; % I-vector calculated by FEKO
Const.FEKOrhsfilename          = 'bowtie-array_MoM.rhs'; % V-vector calculated by FEKO
Const.FEKOoutfilename          = 'bowtie-array_MoM.out';

% --------------------------------------------------------------------------------------------------
% Define output files for transferring expansion coefficients back to FEKO data
% --------------------------------------------------------------------------------------------------
Const.SUNEMdgfmstrfilename     = '';

% --------------------------------------------------------------------------------------------------
% Define additional program flow constants
% --------------------------------------------------------------------------------------------------
% TO-DO: Setup some documentation for this - also assign default values if they are not
% defined explicitely here.
%Const.no_mutual_coupling_array = false; % Deactivate coupling between domains.

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

% Calculate now an interpolated version of the impedance matrix over a
% frequency range. We do not implement the MoM, but rather use FEKO's
% solutions at some points to represent the calculated (i.e. exact values)
zMatricesINTERP = calcInterpolatedZmatrices(Const, Solver_setup, zMatricesFEKO);

% --------------------------------------------------------------------------------------------------
% Run the EM solver
% --------------------------------------------------------------------------------------------------
[SolutionFEKO] = runEMsolvers(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO);

[SolutionINTERP] = runEMsolvers(Const, Solver_setup, zMatricesINTERP, yVectorsFEKO, xVectorsFEKO);
