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

     %zMatricesINTERP;
      frequency = Solver_setup.frequencies.samples;
      freqStart = frequency(1);
      freqEnd = frequency(100);
      stepSize = frequency(2) - frequency(1);
      numFreq = Solver_setup.frequencies.freq_num;
      fstep = 2;
      RWGmBasis = Const.numMoMbasis;
      RWGnBasis = Const.numMoMbasis;
      FrequencySamples = [];
      zMatricesINTERPreal = zeros(510, 50, 510); 
      zMatricesINTERPimaginary = zeros(510, 50, 510);
      inerpolatedreal = [];
      interpolatedImaginary = [];
      Interpolate_Zmn = [];
      interpolatedreal = [];
      interpolatedImaginary = []; 
            
%  


%  for freq = 2:fstep:numFreq
  FrequencySamples = [FrequencySamples; frequency(2:fstep:numFreq)];
   for m = 1:RWGmBasis 
    for n = 1:RWGnBasis
         col = 0;
      if m ~= n
        for freq = 2:fstep:numFreq
            col = col+1;
         %FrequencySamples = [FrequencySamples; frequency(freq)];
         %FrequencysampleUnique(n) = Frequencysamples + linspace(0, 1, length(FrequencySamples))*10E-3
         zMatriceCalculatedReal = real(zMatricesINTERP(m,n,freq));
         zMatricesINTERPreal(n,col,m) = zMatriceCalculatedReal;
         zMatriceCalculatedImaginary = imag(zMatricesINTERP(m,n,freq));
         zMatricesINTERPimaginary(n,col,m) = zMatriceCalculatedImaginary;

       end 
     end
   % row = row+1;
   end
  end

for m = 1:RWGmBasis
    for n = 1:RWGnBasis
        if m ~= n 
      % FrequencySamples = [FrequencySamples; frequency(2:fstep:numFreq)];
      %Apply interpolation
      fq = (freqStart:stepSize:freqEnd);
      a = zMatricesINTERPreal(n,:,m);
      xq = interp1(FrequencySamples,a,fq,"linear");
      interpolatedreal = [interpolatedreal; xq];
      %xq = [];
      b = zMatricesINTERPimaginary(n,:,m);
      yq = interp1(FrequencySamples,b,fq,"linear");
      interpolatedImaginary = [interpolatedImaginary; yq];
      %yq = [];
        end
    end
end
 m = 1;




% 
% 
%    for f = 1:99
%          %store the elements in Zmn (2D) in a struct
%          Zmn = reshape((zInterpolatedsamples(:,f)),[RWGmBasis,RWGnBasis]);  %interpolated data
%          Zmnlist = [Zmnlist; Zmn];
%    end
% 
% 
% %storing the data in 3-D (2x2 matrix)xfrequency
% row = 1;    
% for i = 1:99   %1:numFreq
%     for j = 1:RWGmBasis   
%         for k = 1:RWGnBasis
%             InterpolatedValues(j,k,i) = Zmnlist(row,k);
% 
%         end 
%         row = row+1;
%     end
% end


