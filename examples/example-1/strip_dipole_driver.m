% Author: Danie Ludick (dludick@sun.ac.za)
% Project: Strip dipole antenna simulation
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
Const = sunem_initialise('dipoles',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver              = true;

% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'strip_dipole.mat';
Const.FEKOstrfilename          = 'strip_dipole.str';
Const.FEKOrhsfilename          = 'strip_dipole.rhs';
Const.FEKOoutfilename          = 'strip_dipole.out';
%Const.FEKOefefilename          = 'strip_dipole.efe';
%Const.FEKOffefilename          = 'strip_dipole.ffe';

% The Following file is used to port solutions to FEKO 
% (for post-processing in POSTFEKO).
% TO-DO: [DL] Add this.
% Const.output_strfilename    = '';
% Const.writeFEKOstrfile = [0 0 0 0];

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

% --------------------------------------------------------------------------------------------------
% Postprocess the results, e.g. calculate the Electric field
% --------------------------------------------------------------------------------------------------
ax1 = nexttile;
yvalues = log10(abs(zMatrices.values(1,1,1))); % build 3D array of all of individuals to manipulate as one
yvalues=reshape(permute(yvalues,[1]),1,[]);  % rearrange by plane first, row & column and put in columns
plot(1,yvalues)                                    % plot each column against the y vector

yvalues = log10(abs(zMatrices.values(1,10,1))); 
yvalues=reshape(permute(yvalues,[1]),1,[]);
hold on;
plot(1,yvalues);

 
yvalues = log10(abs(zMatrices.values(1,20,1))); 
yvalues=reshape(permute(yvalues,[1]),1,[]);
hold on;
plot(1,yvalues);

legend('m,n = 1,1','m,n = 1,10','m,n = 1,20');
title(ax1,'magnitude plots');
hold off

