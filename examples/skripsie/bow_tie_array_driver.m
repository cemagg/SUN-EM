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
%[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);

%Solution.mom has all the solver settings

% 
%     f = figure;
%     ax1 = axes(f);
%     hold(ax1,'on');
%     
%     f = figure;
%     ax2 = axes(f);
%     hold(ax2,'on');
    %f = figure;
    %ax3 = axes(f);
    %hold(ax3,'on');

    Zmnn = [];
    Zmnlist = [];
 %   Zmn = reshape(Zmn,[10,10]);
    
% for freq = 1:max 
 for  m = 1:100
  for  n= 1:100
   % if m~= n 

      frequency = Solver_setup.frequencies.samples;
      freqNum = length(Solver_setup.frequencies.samples);
      fstep = 2; %physconst('LightSpeed')/(2*maxRmn)/2;    %18,84MHz
      stepSize = frequency(2) - frequency(1);   %interval between adjacent 'selected' frequencies
      freqStart = frequency(1);
      freqEnd = frequency(Solution.mom.numSols);

      max = Solution.mom.numSols;
      maxRmn = 3.978866829890139;
      newFrequency = frequency(1:fstep:max);   %fewer frequncy samples points
      NewnumSols = length(newFrequency);
      lambda = physconst('LightSpeed')./newFrequency;

      edge_m_X = Solver_setup.rwg_basis_functions_shared_edge_centre(m,1);
      edge_m_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(m,2);

      edge_n_X = Solver_setup.rwg_basis_functions_shared_edge_centre(n,1);
      edge_n_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(n,2);

      Rmn = sqrt((edge_m_X - edge_n_X)^2 + (edge_m_Y - edge_n_Y)^2);
       %if Rmn >= (0.5*lambda)
    
          matrix_Z = zMatrices.values(m,n,1:fstep:max);   % build 3D array of all of individuals to manipulate as one
          matrix_Z = reshape(permute(matrix_Z,[5,4,3,2,1]),NewnumSols,[]); % reshape vector to matrix
          real_z1 = real(matrix_Z);
          imag_z1 = imag(matrix_Z);
    
          % Improved Zmn Solution.mom.numSols
          new_matrixZ = matrix_Z./exp(-1i*((2*pi)./lambda')*Rmn);
    
          new_real1 = real(new_matrixZ);
          new_imag1 = imag(new_matrixZ);
        %  hold on;
        
          % before interpolating the real points
          %plot(ax1,newFrequency,new_real1,'-x');
      
          %Apply interpolation
          fq = (freqStart:stepSize:freqEnd);        
          vq = interp1(newFrequency,new_real1,fq,"spline");
          vr = interp1(newFrequency,new_imag1,fq,"spline");
          %plot(ax1,newFrequency,new_real1,fq,vq,'-b');   
    
%           xlabel(ax1,'FREQUENCY');
%           ylabel(ax1,'RESISTANCE (OHM)');
%           title(ax1,'Real plot');
         
          %before interpolating the imaginary points
          %plot(ax2,newFrequency,new_imag1,'-x');
    
          %Apply interpolation
          %plot(ax2,newFrequency,new_imag1,fq,vr,'-b');  %(Spline interp1)
          
%           xlabel(ax2,'FREQUENCY');
%           ylabel(ax2,'REACTANCE (OHMS)');
%           title(ax2,'Imaginary plot');
    
          %plot(ax3,new_real1,new_imag1,'-x');
          %legends{m,n} = sprintf('m,n = %d,%d', m,n);
         
          %legend( ax1,legends );
          %legend( ax2,legends );
          %hold(ax1,'off');
          %hold(ax2,'off');
    
         %Find error norm percentage between Zinterp and original 
         Zinterp1 = reshape((vq(1:fstep:max) + 1i*(vr(1:fstep:max))),[],1); %reshape vector to matrix
         Zinterp2 = Zinterp1.*exp(-1i*((2*pi)./lambda')*Rmn);             %normalise
         errorNormPercentage = (norm(matrix_Z - Zinterp1)/(norm(new_matrixZ)))* 100;
       

         % now add field called zInterp
         for l = m
          for k = n
             Zmnn = [Zmnn; (vq + 1i*(vr))];
             [zMatrices(:).zInterp] = Zmnn;
          end
         end 

   %end
  %end
 end
 end
   
      for f = 1:99
         
          Zmn = reshape((Zmnn(:,f)),[100,100]);
         Zmnlist = [Zmnlist; Zmn];

         [zMatrices(:).zInterpValues] = Zmn;
         [zMatrices(:).zInterpValues] = Zmnlist;

       end



row = 1;
for i = 1:99

    for j = 1:100
        for k = 1:100
            InterpolatedValues(j,k,i) = Zmnlist(row,k);
            [zMatrices(:).zInterpValues] = InterpolatedValues;

        end 
        row = row+1;
    end
end
[Solution] = runEMsolvers(Const, Solver_setup, zMatrices, yVectors, xVectors);


