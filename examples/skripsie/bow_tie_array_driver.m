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

   m = 1;
   n= 500;
    
      frequency = Solver_setup.frequencies.samples;
      fstep = 5;
      max = Solution.mom.numSols;
      
      newFrequency = frequency(1:fstep:max);
      NewnumSols = length(newFrequency);

      matrix_Z = zMatrices.values(m,n,1:fstep:max);   % build 3D array of all of individuals to manipulate as one
      matrix_Z = reshape(permute(matrix_Z,[5,4,3,2,1]),NewnumSols,[]); % rearrange by plane first, row & column and put in columns
      real_z1 = real(matrix_Z);
      imag_z1 = imag(matrix_Z);
      lambda = physconst('LightSpeed')./newFrequency;
      

      % Improved Zmn Solution.mom.numSols
      edge_m_X = Solver_setup.rwg_basis_functions_shared_edge_centre(m,1);
      edge_m_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(m,2);

      edge_n_X = Solver_setup.rwg_basis_functions_shared_edge_centre(n,1);
      edge_n_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(n,2);
      Rmn = sqrt((edge_m_X - edge_n_X)^2 + (edge_m_Y - edge_n_Y)^2);
    
      new_matrixZ = matrix_Z./exp(-1i*((2*pi)./lambda')*Rmn);
      new_real1 = real(new_matrixZ);
      new_imag1 = imag(new_matrixZ);
      hold on;
    
      % before interpolating the real points
      plot(ax1,newFrequency,new_real1,'-x');
  
      %Apply interpolation
      fq = (frequency(1):12627000:1350270000);            %step size of 200
      vq = interp1(newFrequency,new_real1,fq,"spline");
      vr = interp1(newFrequency,new_imag1,fq,"spline");
      plot(ax1,newFrequency,new_real1,fq,vq,'-b');   

      xlabel(ax1,'FREQUENCY');
      ylabel(ax1,'RESISTANCE (OHM)');
      title(ax1,'Real plot');
     
      %before interpolating the imaginary points
      plot(ax2,newFrequency,new_imag1,'-x');

      %Apply interpolation
      plot(ax2,newFrequency,new_imag1,fq,vr,'-b');  %(Spline interp1)
      
      xlabel(ax2,'FREQUENCY');
      ylabel(ax2,'REACTANCE (OHMS)');
      title(ax2,'Imaginary plot');
      legends{m,n} = sprintf('m,n = %d,%d', m,n);
     
      legend( ax1,legends );
      legend( ax2,legends );
      hold(ax1,'off');
      hold(ax2,'off');


      Zinterp = vq(1:2:max) + 1i*(vr(1:2:max));