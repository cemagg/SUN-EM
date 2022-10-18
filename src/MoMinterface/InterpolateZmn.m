function [Interpolate_Zmn] = InterpolateZmn(Const, Solver_setup, zMatricesINTERP)
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
   
     % Initialise the return values
     zMatricesINTERPp = [];
      zMatricesINTERPp.name = zMatricesINTERPp;
%      numSols = 1;
%      Interpolate_Zmn.numSols = numSols;
     
      FrequencySamples = [];
      zMatricesINTERPreal = zeros(510, 50, 510);  
      zMatricesINTERPimaginary = zeros(510, 50, 510);
      %interpolatedreal = zeros(510, 100);
      %interpolatedImaginary = zeros(510, 100);
      interpolatedreal = [];
      interpolatedImaginary = [];
     Interpolate_Zmn = [];




      frequency = Solver_setup.frequencies.samples;
      freqStart = frequency(1);
      freqEnd = frequency(100);
      stepSize = frequency(2) - frequency(1);
      numFreq = Solver_setup.frequencies.freq_num;
      fstep = 2;
      RWGmBasis = Const.numMoMbasis;
      RWGnBasis = Const.numMoMbasis;
            
%  

%Extract the real and imaginary values preparing for interpolation,
  FrequencySamples = [FrequencySamples; frequency(2:fstep:numFreq)];
   for m = 1:RWGmBasis 
    for n = 1:RWGnBasis
         col = 0;
      if m ~= n
        for freq = 2:fstep:numFreq
            %tic
            col = col+1;
         zMatriceCalculatedReal = real(zMatricesINTERP(m,n,freq)); 
         zMatricesINTERPreal(n,col,m) = zMatriceCalculatedReal;
         zMatriceCalculatedImaginary = imag(zMatricesINTERP(m,n,freq));
         zMatricesINTERPimaginary(n,col,m) = zMatriceCalculatedImaginary;
       %  zMatricesINTERPp.setupTime(freq) = toc;
       end 
      end
      
   % row = row+1;
   end
   end

     %Apply interpolation
 for m = 1:RWGmBasis
    for n = 1:RWGnBasis
       % tic
       % tic
        if m ~= n 
      fq = (freqStart:stepSize:freqEnd);
      a = zMatricesINTERPreal(n,:,m);
      xq = interp1(FrequencySamples,a,fq,"spline");
      interpolatedreal = [interpolatedreal; xq];
      b = zMatricesINTERPimaginary(n,:,m);
      yq = interp1(FrequencySamples,b,fq,"spline");
      interpolatedImaginary = [interpolatedImaginary; yq];
      %zMatricesINTERPp.intepolationTime(n) = toc;
        end
    end
end

     %Interpolatd elements (510x510)-(Diagonals)= 259590 X 99 
     Interpolate_Zmn.values = interpolatedreal(1:fstep:numFreq) + 1i*(interpolatedImaginary(1:fstep:numFreq));



