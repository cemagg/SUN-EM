% Author: Danie Ludick (dludick@sun.ac.za)
% Project: Bowtie array simulation using MoM and DGFM
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
Const = sunem_initialise('bow_tie_array',false);

% --------------------------------------------------------------------------------------------------
% Program flow settings
% --------------------------------------------------------------------------------------------------

% Choose the solvers that will be executed
Const.runMoMsolver          = true;
Const.runDGFMsolver         = false;



% --------------------------------------------------------------------------------------------------
% Define input files for extracting FEKO data
% --------------------------------------------------------------------------------------------------
Const.FEKOmatfilename          = 'bow_tie_array_MoM.mat';
Const.FEKOstrfilename          = 'bow_tie_array_MoM.str';
Const.FEKOrhsfilename          = 'bow_tie_array_MoM.rhs';
Const.FEKOoutfilename          = 'bow_tie_array_MoM.out';


%Const.FEKOstrfilename          = 'bow_tie_array_DGFM.str';
%Const.FEKOrhsfilename          = 'bow_tie_array_DGFM.rhs';
%Const.FEKOoutfilename          = 'bow_tie_array_DGFM.out';



% --------------------------------------------------------------------------------------------------
% Define output files for transferring expansion coefficients back to FEKO data
% --------------------------------------------------------------------------------------------------
Const.SUNEMmomstrfilename      =  ''; %'sunem_mom_bow_tie_array.str';
Const.SUNEMdgfmstrfilename     =  ''; %'sunem_dgfm_bow_tie_array.str';

% --------------------------------------------------------------------------------------------------
% Define additional program flow constants
% --------------------------------------------------------------------------------------------------
% TO-DO: Setup some documentation for this
% Defined explicitely
%Const.no_mutual_coupling_array = false; % Deactivate coupling between domains.

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
% (Note: We either pass our own (internal) matrices, or those read from FEKO). For this particular
% array configuration, we are not yet supporting radiating elements. But as we are consistent with the
% FEKO mesh, we can just use the FEKO RHS vector.
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);

%Solution.mom has all the solver settings


%Top plot
ax1 = nexttile;
frequency = Solver_setup.frequencies.samples;
matrix_Z = zMatrices.values(1,200,1:85);   % build 3D array of all of individuals to manipulate as one
matrix_Z = reshape(permute(matrix_Z,[5,4,3,2,1]),85,[]); % rearrange by plane first, row & column and put in columns
real_z1 = real(matrix_Z);
imag_z1 = imag(matrix_Z);

%
edge_m_X = Solver_setup.rwg_basis_functions_shared_edge_centre(1,1);
edge_m_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(1,2);

edge_n_X = Solver_setup.rwg_basis_functions_shared_edge_centre(200,1);
edge_n_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(200,2);

Rmn = sqrt((edge_m_X - edge_n_X)^2 + (edge_m_Y - edge_n_Y)^2);


new_matrixZ = matrix_Z/exp(-1i*2*pi*Rmn);
new_real1 = real(new_matrixZ);
new_imag1 = imag(new_matrixZ);
hold on;

plot(frequency,new_real1,'-x');
plot(frequency,real_z1,'-x');

%Apply interpolation
%xq = (100131000:200:1350270000);
%Interp1 = spline(frequency,new_real1,xq);
%Interp2 = spline(frequency,new_imag1,xq);
%plot(frequency,new_real1,'-x',xq,Interp1);
%plot(frequency,real_z1,'-x',xq,Interp1);



xlabel('FREQUENCY');
ylabel('RESISTANCE (OHM)');
legend('Improved','old');
title(ax1,'Real plot');
hold off;

ax2 = nexttile;
hold on;

plot(frequency,new_imag1,'-x');
plot(frequency,imag_z1,'-x');

%Apply interpolation
%plot(frequency,new_imag1,'-x',xq,Interp2);
%plot(frequency,imag_z1,'-x',xq,Interp2);

xlabel('FREQUENCY');
ylabel('REACTANCE (OHMS)');
legend('Improved', 'old');
title(ax2,'Imaginary plot');
hold off; 