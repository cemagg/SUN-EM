function [Interpolated_Data] = InterpolateZmn(Const, Solver_setup, zMatricesFEKOcalulated, zMatricesFEKO)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details        
    %       yVectors
    %           The Yrhs-vector data
    %
    %   Output Arguments:
    %       Interpolated_Zmn
    %           Structs containing Zmn calculated after interpolation and timing data
    %
    %   Description:
    %   interpolate Zmn

    % zMatricesInterpolated
    %


    narginchk(4,4);
   
     % Initialise the return values
      Interpolated_Data.name = 'Interpolated_Data';
      Interpolated_Data.Interpolate_Zmn = [];
      Interpolated_Data.interpTime = [];
      Interpolated_Data.numFreq = [];
      Interpolated_Data.errorNormPercentage = [];
      RWGmBasis = Solver_setup.num_mom_basis_functions;
      RWGnBasis = Solver_setup.num_mom_basis_functions;
      Interpolated_Data.Interpolate_ZmnValues = zeros(size(zMatricesFEKO.values));
      Interpolated_Data.errorNormPercentages = [];
      z = [];
      y = [];
      error = [];
     

      frequency = Solver_setup.frequencies.samples;
      stepSize = frequency(2) - frequency(1);
      freqStart = frequency(1);
      freqEnd = frequency(3) + stepSize;
      numFreq = Solver_setup.frequencies.freq_num;
      fstep = 2;


      FrequencySamples = [];
      zMatricesINTERPreal = zeros(RWGmBasis, RWGmBasis);  
      zMatricesINTERPimaginary = zeros(RWGmBasis, RWGmBasis);
      Interpolate_Zmn = zeros(size(zMatricesFEKO.values)); 
      interpTime = [];
      sum = 0;
      sums = [];

%Extract the real and imaginary values preparing for interpolation,
  %FrequencySamples = [FrequencySamples; frequency(2:fstep:numFreq)];
  for freq = 1:fstep:numFreq
   for m = 1:RWGmBasis 
    for n = 1:RWGnBasis
           zMatricesINTERPreal(n,m,freq) = real(zMatricesFEKOcalulated(n,m,freq));
           zMatricesINTERPimaginary(n,m,freq) = imag(zMatricesFEKOcalulated(n,m,freq));
    end
   end
  end


for m = 1:RWGmBasis
   % tic
   for n = 1:RWGmBasis 
     if m ~= n
      for freq = 1:3  
         z = [z, zMatricesINTERPreal(n,m,freq)];
         y = [y, zMatricesINTERPimaginary(n,m,freq)];
       end
         NewfrequencySamples = frequency(1:fstep:numFreq);
         zs = z == 0;
         z(zs) = [];
         xs = y == 0;
         y(xs) = [];
         if size(NewfrequencySamples) == numel(z)
           fq = (freqStart:stepSize:freqEnd);
           xq = interp1(NewfrequencySamples,z,fq,"spline");
           yq = interp1(NewfrequencySamples,y,fq,"spline");
             for i = 1:numFreq
                 threeD(n,m,i) = xq(1,i);
                 threeDd(n,m,i) = yq(1,i);
                 Interpolate_Zmn(n,m,i) = threeD(n,m,i) + 1i*threeDd(n,m,i);
             end
         end
       z = [];
       y = [];
     else 
       Interpolate_Zmn(n,m,:) = zMatricesFEKOcalulated(n,m,:);
     end
   end
end




 %nomramlise the interpolated Zmn
 for freq = 1:3
    for m = 1:RWGmBasis 
        for n = 1:RWGnBasis
            if m ~= n
                FrequencySample = frequency(freq);
                lambda = physconst('LightSpeed')/FrequencySample;

                edge_m_X = Solver_setup.rwg_basis_functions_shared_edge_centre(m,1);
                edge_m_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(m,2);

                edge_n_X = Solver_setup.rwg_basis_functions_shared_edge_centre(n,1);
                edge_n_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(n,2);

                Rmn = sqrt((edge_m_X - edge_n_X)^2 + (edge_m_Y - edge_n_Y)^2);

               Interpolate_Zmn(n,m,freq) = Interpolate_Zmn(n,m,freq)*exp(-1i*((2*pi)/lambda')*Rmn);
            end
            Interpolated_Data.Interpolate_ZmnValues(n,m,:) = Interpolate_Zmn(n,m,:);
        end 
    end
 end

 
for freq = 1:3
    refMatrix = complex(zeros(size(zMatricesFEKO.values)));
    refMatrix = refMatrix + zMatricesFEKO.values(:,freq);

   %errorNormPercentage = calculateMatrixErrorNormPercentage(refMatrix, refMatrix);
   % fprintf(1,"Frobenius Error Norm Percentage = %f percent\n",errorNormPercentage);

    newMatrix = refMatrix +  Interpolated_Data.Interpolate_ZmnValues(:,freq)./10000;

    errorNormPercentage = calculateMatrixErrorNormPercentage(refMatrix, newMatrix);
    fprintf(1,"Frobenius Error Norm Percentage = %f percent\n",errorNormPercentage);
end


%  
% Example to show how to calculate the frobenius norm associated with
% 2 x complex matrices, where 1 is the reference solution and the other
%  one that we want to figure out what the accuracy is​
% Generate some reference data.
% for freq = 1:99
%     refMatrix = complex(zeros(size(zMatrices.values)));
%     refMatrix = refMatrix + zMatrices.values(:,freq);​
%     errorNormPercentage = calculateMatrixErrorNormPercentage(refMatrix, refMatrix);
%     fprintf(1,"Frobenius Error Norm Percentage = %f percent\n",errorNormPercentage);
% 
%    newMatrix = refMatrix +  Interpolated_Data.Interpolate_ZmnValues(:,freq)./1000;
% 
%    errorNormPercentage = calculateMatrixErrorNormPercentage(refMatrix, newMatrix);
%    fprintf(1,"Frobenius Error Norm Percentage = %f percent\n",errorNormPercentage);
% end
 q = 1;



% 
% Now generate another matrix (e.g. with another technique used to
% calculate the matrix elements). Here we just use the same values as in
% refMatrix, but we add some small random numbers (you can adjust the
% accuracy by adjusting the division factor -> e.g. 1000 will result in a
% very small Frob norm (i.e. the matrix elements are essentailly the same -
% verify this by looking at the elements).
%  newMatrix = refMatrix + rand(3)./100 + 1i.*rand(3)./100;
% % ​
% % Calculate the error norm percentage of the new matrix by comparing it
% % against the reference data.
% errorNormPercentage = calculateMatrixErrorNormPercentage(refMatrix, newMatrix);
% fprintf(1,"Frobenius Error Norm Percentage = %f percent\n",errorNormPercentage);
% % 









% sum = 0;
% for freq=1:99
%        Interpolated_Data.errorNormPercentage = (norm(zMatrices.values(:,freq) - Interpolate_Zmn(:,freq))/norm(zMatrices.values(:,freq)))* 100;
%        error = [error; ((norm(zMatrices.values(:,freq) - Interpolate_Zmn(:,freq))/norm(zMatrices.values(:,freq)))* 100)];
%        sum = sum + ((norm(zMatrices.values(:,freq) - Interpolate_Zmn(:,freq))/norm(zMatrices.values(:,freq)))* 100);
%        Interpolated_Data.errorNormPercentages(n,m) = (zMatrices.values(:,freq) - Interpolate_Zmn(:,freq))/(zMatrices.values(:,freq))* 100;
%  sums = [sums; sum];
%  sum = 0;
% end





