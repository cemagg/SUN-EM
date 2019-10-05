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
Const = sunem_initialise('Close_array_active',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver          = true;
Const.runCBFMsolver         = true;

% --------------------------------------------------------------------------------------------------
% Define input array files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'Close_array_active.mat';
Const.FEKOstrfilename          = 'Close_array_active.str';
Const.FEKOrhsfilename          = 'Close_array_active.rhs';
Const.FEKOoutfilename          = 'Close_array_active.out';
Const.FEKOefefilename          = 'Close_array_active.efe';
Const.FEKOffefilename          = 'Close_array_active.ffe';

% --------------------------------------------------------------------------------------------------
% Define output files for transferring expansion coefficients back to FEKO data
% --------------------------------------------------------------------------------------------------
Const.SUNEMmomstrfilename      =  '';%'sunem_mom_bow_tie_array.str';

% --------------------------------------------------------------------------------------------------
% Define additional program flow constants
% --------------------------------------------------------------------------------------------------
Const.no_mutual_coupling_array = false; % Deactivate coupling between domains.
Const.calcSecMBFs = true;               % For MBF based solvers
Const.useMBFreduction = true;           % SVD applied after the MBFs are generated to retain an orthonormal set
Const.MBFthreshold = 1000;              % Threshold used for the SVD reduction of the MBFs

% --------------------------------------------------------------------------------------------------
% Read the MoM full array matrix equation from the file
% --------------------------------------------------------------------------------------------------
[Const, zArrayMatricesFEKO, yArrayVectorsFEKO, xArrayVectorsFEKO] = extractFEKOMoMmatrixEq(Const);


% --------------------------------------------------------------------------------------------------
% Parse the setup files to extract the frequency sweep, the geometry and basis function setup 
% --------------------------------------------------------------------------------------------------
% TO-DO: At a later stage we can also add other meshing / geometry
% preprocessxing, e.g. Gmsh or GiD. For now the solver setup is read from FEKO.
[Const, Solver_setup] = parseFEKOoutfile(Const, yArrayVectorsFEKO);
% -------------------------------------------------------------------------------------------------
% Change the variables below to change what is included in the MBF
% generation
% -------------------------------------------------------------------------------------------------
SourceMode = 0;
DGFM = false;

% --------------------------------------------------------------------------------------------------
% Run the CMAMBFs generator
% --------------------------------------------------------------------------------------------------
% This for loop adjusts the amount of modes used as MBFs the timing results
% have been taken out temporarily, only the error norm is displayed for
% each loop.
for i=5:27
    [CMAMBFs] = runCMA_MBFgenerator(Const, Solver_setup, zArrayMatricesFEKO, yArrayVectorsFEKO, DGFM,SourceMode, i);
    [CBFMCMA] = runCMACBFM(Const, Solver_setup, zArrayMatricesFEKO, yArrayVectorsFEKO, xArrayVectorsFEKO, CMAMBFs);
end

