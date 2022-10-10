f = figure;
ax1 = axes(f);
hold(ax1,'on');

f = figure;
ax2 = axes(f);
hold(ax2,'on');

for m = 1
 for n = 200
 % if m == n
      
     % ax1 = nexttile;
      frequency = Solver_setup.frequencies.samples;
      %frequency_step_size = frequency ;
      newFrequency = Solver_setup.frequencies.samples(1:2:100);
      NewnumSols = length(Solver_setup.frequencies.samples(1:2:100));
    
      matrix_Z = zMatrices.values(m,n,1:NewnumSols);   % build 3D array of all of individuals to manipulate as one
      matrix_Z = reshape(permute(matrix_Z,[5,4,3,2,1]),NewnumSols,[]); % rearrange by plane first, row & column and put in columns
      real_z1 = real(matrix_Z);
      imag_z1 = imag(matrix_Z);
      lambda = physconst('LightSpeed')./newFrequency;


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
      vq = interp1(newFrequency,new_real1,fq,"spline");
      vr = interp1(newFrequency,new_imag1,fq,"spline");
      plot(ax1,newFrequency,new_real1,fq,vq,'-');
    
      xlabel(ax1,'FREQUENCY');
      ylabel(ax1,'RESISTANCE (OHM)');
      title(ax1,'Real plot');
    

     % before interpolating the imaginary points
      %ax2 = nexttile;
      %plot(ax2,frequency,imag_z1,'-xb');
      plot(ax2,newFrequency,new_imag1,'-xb');


      %Apply interpolation
      plot(ax2,newFrequency,new_imag1,fq,vr,'-');
      %plot(frequency,imag_z1,'-x',fq,Interp2);
     
      xlabel(ax2,'FREQUENCY');
      ylabel(ax2,'REACTANCE (OHMS)');
      title(ax2,'Imaginary plot');
      %legends{m,n} = sprintf('m,n = %d,%d', m,n);

  %end
 end
end
%legend(ax2,'Improved sample points');
 %legend( ax1,legends );
 %legend( ax2,legends );

hold(ax1,'off');

hold(ax2,'off');