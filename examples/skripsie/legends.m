%Top plot
ax1 = nexttile;
frequency = Solver_setup.frequencies.samples;
matrix_Z = zMatrices.values(1,1,1:Solution.mom.numSols);   % build 3D array of all of individuals to manipulate as one
matrix_Z = reshape(permute(matrix_Z,[5,4,3,2,1]),Solution.mom.numSols,[]); % rearrange by plane first, row & column and put in columns
real_z1 = real(matrix_Z);
imag_z1 = imag(matrix_Z);
lambda = physconst('LightSpeed')./frequency;

% 
edge_m_X = Solver_setup.rwg_basis_functions_shared_edge_centre(1,1);
edge_m_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(1,2);

edge_n_X = Solver_setup.rwg_basis_functions_shared_edge_centre(1,1);
edge_n_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(1,2);

Rmn = sqrt((edge_m_X - edge_n_X)^2 + (edge_m_Y - edge_n_Y)^2);


new_matrixZ = matrix_Z./exp(-1i*((2*pi)./lambda')*Rmn);
new_real1 = real(new_matrixZ);
new_imag1 = imag(new_matrixZ);
hold on;

%plot(frequency,real_z1,'-xb');
plot(frequency,new_real1,'-xr');


%Apply interpolation
%fq = (100131000:200:1350270000);            %step size of 200
%Interp1 = spline(frequency,new_real1,fq);
%Interp2 = spline(frequency,new_imag1,fq);
%plot(frequency,new_real1,'xr',fq,Interp1,'-b');
%plot(frequency,real_z1,fq,Interp2,'-xb');
%legend('Original','Improved sample points','spline');



xlabel('FREQUENCY');
ylabel('RESISTANCE (OHM)');
legend('Improved sample points');
title(ax1,'Real plot');
hold off;

ax2 = nexttile;
hold on;

%plot(frequency,imag_z1,'-xb');
plot(frequency,new_imag1,'-xr');


%Apply interpolation
%plot(frequency,new_imag1,'xr',fq,Interp2,'-r');
%plot(frequency,imag_z1,'-x',fq,Interp2);

xlabel('FREQUENCY');
ylabel('REACTANCE (OHMS)');
legend('Improved sample points');
title(ax2,'Imaginary plot');
hold off; 





for i = 1:1:zMatrices.mBasis
  for j = 1:1:zMatrices.nBasis