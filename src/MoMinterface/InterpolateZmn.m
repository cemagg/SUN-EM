function [Interpolated_Data] = InterpolateZmn(Const, Solver_setup, zMatricesINTERP, zMatrices)
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
      RWGmBasis = Solver_setup.num_mom_basis_functions; %Solver_setup.num_mom_basis_functions;
      RWGnBasis = Solver_setup.num_mom_basis_functions; %Solver_setup.num_mom_basis_functions;
      Interpolated_Data.Interpolate_ZmnValues = zeros(RWGmBasis, RWGnBasis, 100);
      z = [];
      y = [];
     

      frequency = Solver_setup.frequencies.samples;
      stepSize = frequency(2) - frequency(1);
      freqStart = frequency(1);
      freqEnd = frequency(100) + stepSize;
      numFreq = Solver_setup.frequencies.freq_num;
      fstep = 2;

%       RWGmBasis = 510; %Const.numMoMbasis;
%       RWGnBasis = 510; %Const.numMoMbasis;

      FrequencySamples = [];
      zMatricesINTERPreal = zeros(RWGmBasis, RWGmBasis);  
      zMatricesINTERPimaginary = zeros(RWGmBasis, RWGmBasis);
      Interpolate_Zmn = zeros(RWGmBasis, RWGmBasis, 100); 
      interpTime = [];

%Extract the real and imaginary values preparing for interpolation,
  %FrequencySamples = [FrequencySamples; frequency(2:fstep:numFreq)];
  for freq = 1:fstep:numFreq
   for m = 1:RWGmBasis 
    for n = 1:RWGnBasis
         col = 0;
      %if m ~= n
        %for freq = 2:fstep:numFreq
           col = col+1;
           zMatriceCalculatedReal = real(zMatricesINTERP(n,m,freq)); 
           zMatricesINTERPreal(n,m,freq) = zMatriceCalculatedReal;
           zMatriceCalculatedImaginary = imag(zMatricesINTERP(n,m,freq));
           zMatricesINTERPimaginary(n,m,freq) = zMatriceCalculatedImaginary;
      % end 
     % end
    end
   end
  end


for m = 1:RWGmBasis
   for n = 1:RWGmBasis 
    for freq = 1:99
         B = zMatricesINTERPreal(n,m,freq);
         z = [z, B];
         C = zMatricesINTERPimaginary(n,m,freq);
         y = [y, C];
    end
         NewfrequencySamples = frequency(1:fstep:numFreq);
         zs = z == 0;
         z(zs) = [];
         xs = y == 0;
         y(xs) = [];
         fq = (freqStart:stepSize:freqEnd);
         xq = interp1(NewfrequencySamples,z,fq,"spline");
         yq = interp1(NewfrequencySamples,y,fq,"spline");
         for i = 1:numFreq
        threeD(n,m,i) = xq(1,i);
        threeDd(n,m,i) = yq(1,i);
        Interpolate_Zmn(n,m,i) = threeD(n,m,i) + 1i*threeDd(n,m,i);
         end
       z = [];
       y = [];
  end
end




%  %Apply interpolation and store interpolated data in Zmn
%  for m = 1:RWGmBasis
%     for n = 1:RWGnBasis
%         tic
%         if m ~= n 
%             fq = (freqStart:stepSize:freqEnd);
%             a = zMatricesINTERPreal(n,m);
%             xq = interp1(frequency,a,fq,"spline");
%             b = zMatricesINTERPimaginary(n,m,:);
%             yq = interp1(frequency,b,fq,"spline");
%             Interpolate_Zmn(n,m,:) = xq + 1i*yq;       %add to get Zmn = xq + j*yq
%             Interpolated_Data.interpTime(m) = toc;
%         end
%     end
%  end



 %nomramlise the interpolated Zmn
 for freq = 1:fstep:numFreq
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
               Interpolated_Data.Interpolate_ZmnValues(n,m,:) = Interpolate_Zmn(n,m,:);

           end
        end 
    end
 end
% m = 1;
 

for freq=1:fstep:numFreq
 for m = 1:RWGmBasis
   for n = 1:RWGnBasis
      if m ~=n
       Interpolated_Data.errorNormPercentage(n,m) = (norm(zMatrices.values(n,m,freq) - norm(Interpolate_Zmn(n,m,freq)))/(norm(zMatrices.values(n,m,freq))))* 100; 
      end
   end
  end
end
q = 1;




