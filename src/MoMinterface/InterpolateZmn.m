function [Interpolate_Zmn] = InterpolateZmn(Const, Solver_setup, zMatricesINTERP, zMatrices)
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
      Interpolated_Data = [];
      Interpolated_Data.name = Interpolated_Data;
      Interpolated_Data.Interpolate_Zmn = [];
      Interpolated_Data.interpTime = [];
      Interpolated_Data.numFreq = [];
      Interpolated_Data.errorNormPercentage = [];
     

      frequency = Solver_setup.frequencies.samples;
      stepSize = frequency(2) - frequency(1);
      freqStart = frequency(1);
      freqEnd = frequency(100) + stepSize;
      numFreq = Solver_setup.frequencies.freq_num;
      fstep = 2;
      RWGmBasis = 200; %Const.numMoMbasis;
      RWGnBasis = 200; %Const.numMoMbasis;

      FrequencySamples = [];
      zMatricesINTERPreal = zeros(RWGmBasis, 50, RWGmBasis);  
      zMatricesINTERPimaginary = zeros(RWGmBasis, 50, RWGmBasis);
      Interpolate_Zmn = zeros(RWGmBasis, RWGmBasis, 100);           

%Extract the real and imaginary values preparing for interpolation,
  FrequencySamples = [FrequencySamples; frequency(2:fstep:numFreq)];
   for m = 1:RWGmBasis 
    for n = 1:RWGnBasis
         col = 0;
      if m ~= n
        for freq = 2:fstep:numFreq
            %tic
           col = col+1;
           zMatriceCalculatedReal = real(zMatricesINTERP(n,m,freq)); 
           zMatricesINTERPreal(n,col,m) = zMatriceCalculatedReal;
           zMatriceCalculatedImaginary = imag(zMatricesINTERP(n,m,freq));
           zMatricesINTERPimaginary(n,col,m) = zMatriceCalculatedImaginary;
          %zMatricesINTERPp.setupTime(freq) = toc;
       end 
      end
   end
   end

 %Apply interpolation and store interpolated data in Zmn
 for m = 1:RWGmBasis
    for n = 1:RWGmBasis
       % tic
       % tic
        if m ~= n 
            fq = (freqStart:stepSize:freqEnd);
            a = zMatricesINTERPreal(n,:,m);
            xq = interp1(FrequencySamples,a,fq,"spline");
            b = zMatricesINTERPimaginary(n,:,m);
            yq = interp1(FrequencySamples,b,fq,"spline");
            Interpolate_Zmn(n,m,:) = xq + 1i*yq;       %add to get Zmn = xq + j*yq
            %zMatricesINTERPp.intepolationTime(n) = toc;
        end
    end
 end


 %nomramlise the interpolated Zmn
 for freq = 1:fstep:numFreq
    for m = 1:RWGmBasis  
        for n = 1:RWGmBasis
            if m ~= n
                
                FrequencySample = frequency(freq);
                lambda = physconst('LightSpeed')./FrequencySample;

                edge_m_X = Solver_setup.rwg_basis_functions_shared_edge_centre(m,1);
                edge_m_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(m,2);

                edge_n_X = Solver_setup.rwg_basis_functions_shared_edge_centre(n,1);
                edge_n_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(n,2);

                Rmn = sqrt((edge_m_X - edge_n_X)^2 + (edge_m_Y - edge_n_Y)^2);

                Interpolate_Zmn(n,m,freq) = Interpolate_Zmn(n,m,freq)*exp(-1i*((2*pi)./lambda')*Rmn);
           end
        end 
    end
 end

for freq=1:fstep:numFreq
 for m = 1:RWGmBasis 
   for n = 1:RWGmBasis
      if m ~=n
           errorNormPercentage(n,m) = (norm(zMatrices.values(n,m,freq) - Interpolate_Zmn(n,m,freq))/(norm(zMatrices.values(n,m,freq))))* 100; 
      end
   end
  end
end
q = 1;

%interpolatedreal = zeros(510, 100);
%interpolatedImaginary = zeros(510, 100);



