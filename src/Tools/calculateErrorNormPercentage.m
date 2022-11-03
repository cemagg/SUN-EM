function errorNormPercentage = calculateErrorNormPercentage(refData, data)
    %calculateErrorNormPercentage
    %   Usage:
    %           fc_init(Const, yVectors)
    %
    %   Input Arguments:
    %       ref, data
    %           An N-dimensional vector for which we will calculate the 2-norm error
    %
    %   Output Arguments:
    %       errorNormPercentage
    %           The relative (2-norm) error expressed as a percentage
    %
    %   Description:
    %       Calculates the relative error norm [%], based on the two input vectors:
    %   =======================

    error(nargchk(2,2,nargin));

    if ( length(refData) ~= length(data) )
        error('[calculateErrorNormPercentage] Data-sets not the same size'); 
    end

    differenceSum = 0;
    referenceSum = 0;

    for n = 1:length(refData)
        differenceSum = differenceSum + (abs( data(n) - refData(n) ))^2;
        referenceSum  = referenceSum  + (abs( refData(n) ))^2;
    end

    errorNormPercentage = (sqrt(differenceSum) / sqrt(referenceSum)) * 100 ;