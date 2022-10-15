function [Interpolate_Zmn] = InterpolateZmn(Const, Solver_setup, zMatricesINTERP)
%InterpolateZmn(Const, Solver_setup, yVectors, Solution.dgfm);
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
    %       Interpolate_Zmn
    %           Structs containing Zmn calculated after interpolation and timing data
    %
    %   Description:
    %   interpolate Zmn

    % zMatricesInterpolated
    %


    narginchk(3,3);

      frequency = Solver_setup.frequencies.samples;
      freqStart = frequency(1);
      freqEnd = frequency(100);
      stepSize = frequency(2) - frequency(1);
      numFreq = Solver_setup.frequencies.freq_num;
      fstep = 2;
      RWGmBasis = Const.numMoMbasis;
      RWGnBasis = Const.numMoMbasis;
      FrequencySamples = [];
   
   for m = 1:RWGmBasis   
     for n = 1:RWGnBasis
         if m ~= n
              for freq = 1:fstep:numFreq
              FrequencySample = frequency(freq);
              FrequencySamples = [FrequencySamples; FrequencySample];
              end

          zMatricesINTERP = zMatricesINTERP(m,n,1:2:numFreq);   % build 3D array of all of individuals to manipulate as one
          zMatricesINTERP = reshape(permute(zMatricesINTERP,[3,2,1]),50,[]); % reshape vector to matrix
          real = real(zMatricesINTERP);
          imaginary = imag(zMatricesINTERP);

          fq = (freqStart:stepSize:freqEnd);        
          xq = interp1(FrequencySamples,real,fq,"spline");
          yq = interp1(FrequencySamples,imaginary,fq,"spline");

         Interpolate_Zmn = (xq+1i*yq); %reshape vector to matrix
         Interpolate_Zmn = Interpolate_Zmn.*exp(-1i*((2*pi)./lambda')*Rmn);               %normalise
         errorNormPercentage = (norm(zMatricesINTERP - Interpolate_Zmn)/(norm(zMatricesINTERP)))* 100; 

         error = [error; errorNormPercentage];

      
        end 
        %zMatricesCalculated(m,n,freq) = [zMatricesCalculated(m,n,freq); zMatricesINTERP(m,n,freq)];
    end
  end



