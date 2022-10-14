function [zMatricesINTERP] = calcInterpolatedZmatrices(Const, Solver_setup, zMatricesFEKO)
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
    zMatricesINTERP  = zMatricesFEKO;
   
    % The code below is just an example:
    
    for m = 1:totRWGbasisFunctions
        for n = 1:totRWGbasisFunctions
            
            % Do not do this for self-interaction matrices
            if (m ~= n)
                % Extract the ALL values of 
                zMatricesValuesMN = zMatricesFEKO.values(m,n,1:totFreqSamples);

                zMatricesValuesMNCalculated = zMatricesValuesMN(1:step:freqSample)

                % The following is probably a separate function that performs
                % the spline interpolation .. ?
                zMatricesValuesMNInterpolated = InterpolateZmn(zMatricesValuesMN(1:step:freqSample)

                % Make sure you set your return values again correct
                zMatricesINTERP.values(m,n,1:totFreqSamples) = zMatricesValuesMNInterpolated
            
        end
            
    