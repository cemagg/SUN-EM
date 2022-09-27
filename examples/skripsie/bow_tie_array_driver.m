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

f = figure;
ax1 = axes(f);
hold(ax1,'on');

f = figure;
ax2 = axes(f);
hold(ax2,'on');

for m = 1:6
 for n = 1:6
     if m == n
    
 % ax1 = nexttile;
    frequency = Solver_setup.frequencies.samples;
    matrix_Z = zMatrices.values(m,n,1:Solution.mom.numSols);   % build 3D array of all of individuals to manipulate as one
    matrix_Z = reshape(permute(matrix_Z,[5,4,3,2,1]),Solution.mom.numSols,[]); % rearrange by plane first, row & column and put in columns
    real_z1 = real(matrix_Z);
    imag_z1 = imag(matrix_Z);
    lambda = physconst('LightSpeed')./frequency;

  % 
    edge_m_X = Solver_setup.rwg_basis_functions_shared_edge_centre(m,1);
    edge_m_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(m,2);

    edge_n_X = Solver_setup.rwg_basis_functions_shared_edge_centre(n,1);
      edge_n_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(n,2);
     
      Rmn = sqrt((edge_m_X - edge_n_X)^2 + (edge_m_Y - edge_n_Y)^2);
    
    
      new_matrixZ = matrix_Z./exp(-1i*((2*pi)./lambda')*Rmn);
      new_real1 = real(new_matrixZ);
      new_imag1 = imag(new_matrixZ);
      hold on;
    
      %plot(ax1,frequency,real_z1,'-xb');
      %plot(ax1,frequency,new_real1,'-xr');

    
      %Apply interpolation
      fq = (100131000:200:1350270000);            %step size of 200
      vq = interp1(frequency,new_real1,fq,"makima");
      vr = interp1(frequency,new_imag1,fq,"makima");
      plot(ax1,frequency,new_real1,'xr',fq,vq,'-b');
      %plot(ax1,frequency,real_z1,fq,Interp2,'-xb');
      %legend('Original','Improved sample points','spline');
    
     
    
      xlabel('FREQUENCY');
      ylabel('RESISTANCE (OHM)');
      legend(ax1,'Improved sample points');
      title(ax1,'Real plot');
    
    
      %ax2 = nexttile;
      %plot(ax2,frequency,imag_z1,'-xb');
      %plot(ax2,frequency,new_imag1,'-xb');


      %Apply interpolation
      plot(ax2,frequency,new_imag1,'xr',fq,vr,'-b');
      %plot(frequency,imag_z1,'-x',fq,Interp2);
     
      xlabel('FREQUENCY');
      ylabel('REACTANCE (OHMS)');
      legend(ax2,'Improved sample points');
      title(ax2,'Imaginary plot');
     end
  end
end
legend show
hold(ax1,'off');

hold(ax2,'off');