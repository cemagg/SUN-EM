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
Const = sunem_initialise('bow_tie_array-3x1',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver       = true;
Const.runDGFMsolver      = true;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'bow_tie_array-3x1.mat'; % Z-matrix calculated by FEKO
Const.FEKOstrfilename          = 'bow_tie_array-3x1.str'; % I-vector calculated by FEKO
Const.FEKOrhsfilename          = 'bow_tie_array-3x1.rhs'; % V-vector calculated by FEKO
Const.FEKOoutfilename          = 'bow_tie_array-3x1.out';

% --------------------------------------------------------------------------------------------------
% Define output files for transferring expansion coefficients back to FEKO data
% --------------------------------------------------------------------------------------------------
Const.SUNEMdgfmstrfilename     = '';

% --------------------------------------------------------------------------------------------------
% Define additional program flow constants
% --------------------------------------------------------------------------------------------------
% Do not change these values.
%Const.no_mutual_coupling_array = false; % Deactivate coupling between domains.
Const.useDGFMmethod = 1;                 % 1: Local matrices.
Const.DGFMweightVectorCalcScheme = 0;    % 0: Uniform coefficients / array excitation law. 

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
%zMatricesINTERP = calcInterpolatedZmatrices(Const, Solver_setup, zMatricesFEKO);

% --------------------------------------------------------------------------------------------------
% Run the EM solver
% --------------------------------------------------------------------------------------------------
[SolutionFEKO] = runEMsolvers(Const, Solver_setup, zMatricesFEKO, yVectorsFEKO, xVectorsFEKO);

%[SolutionINTERP] = runEMsolvers(Const, Solver_setup, zMatricesINTERP, yVectorsFEKO, xVectorsFEKO);
