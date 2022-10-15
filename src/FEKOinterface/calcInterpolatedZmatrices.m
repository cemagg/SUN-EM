function [zMatricesINTERP] = calcInterpolatedZmatrices(Const, Solver_setup, zMatrices)
    %runEMsolvers
    %   Usage:
    %       [zMatricesINTERP] = calcInterpolatedZmatrices(Const, Solver_setup, zMatricesFEKO)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containging settings of which solver to run
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details
    %       zMatrices
    %           The Z-matrices data (all the values from e.g. FEKO)
    %   Output Arguments:
    %      zMatricesINTERP
    %           The Z-matrices data (that includes both the calculated
    %           values from e.g. FEKO, as well as interpolated values
    %
    %   Description:
    %       Reads in impedance Z-matrices over a frequency range. For a
    %       givem Z(m,n) element it then extracts a number of calculated
    %       points and determines the rest using interpolation.

    narginchk(3,3);

    % Initialise the return value to be the same as that of FEKO (to get
    % the general structure correctly setup.

  %%%%%%%  zMatricesINTERP  = zMatricesFEKO;

    % The code below is just an example:

   zMatricesINTERP = [];


    % for freq = 1:numFreq 
    RWGmBasis = 10;
    RWGnBasis = 10;
    
% Zmn = zMatrices.values(1,200,1:5);
for  m = 1:RWGmBasis
  for  n= 1:RWGnBasis
        
    if m ~= n
        
          numFreq = Solver_setup.frequencies.freq_num;
          
          zMatricesFEKO = zMatrices.values(m,n,1:numFreq);   % build 3D array of all of individuals to manipulate as one
          zMatricesFEKO = reshape(permute(zMatricesFEKO,[3,2,1]),numFreq,[]); % reshape vector to matrix

          zMatricesINTERP = [zMatricesINTERP; zMatricesFEKO];

         %Calculated zMatrices at selected frequencies
        % fstep = 2;  %50 evenly spaced frequencies

        % zMatricesFEKOcalulated = zMatrices.values(m,n,1:fstep:numFreq);
         
         
         %call the function for interpolation
        % zMatricesInterpolated = InterpolateZmn(zMatricesFEKOcalulated);
    end
   end 
end
