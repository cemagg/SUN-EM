function errorNormPercentage = calculateMatrixErrorNormPercentage(refData, data)
    %calculateErrorNormPercentage
    %   Usage:
    %           errorNormPercentage = calculateMatrixErrorNormPercentage(refData, data)
    %
    %   Input Arguments:
    %       ref, data
    %           MxN-dimensional matrices (one to be used as reference, the other as 
    %           the matrix that we want to check)
    %
    %   Output Arguments:
    %       errorNormPercentage
    %           The relative (frobenius-norm) error expressed as a percentage
    %
    %   Description:
    %       Calculates the relative error norm [%], based on the two input matrices.
    %       We use an expression similar to the 2-norm relative error percentage for vectors:
    %   =======================

    error(nargchk(2,2,nargin));

    if ( size(refData) ~= size(data) )
        error('[calculateMatrixErrorNormPercentage] Data-sets not the same size'); 
    end

    differenceSum = 0;
    referenceSum = 0;

    num_rows = size(refData,1);
    num_cols = size(refData,2);

    for m = 1:num_rows
        for n = 1:num_cols
            differenceSum = differenceSum + (abs( data(m,n) - refData(m,n) ))^2;
            referenceSum  = referenceSum  + (abs( refData(m,n) ))^2;
        end %for n = 1:num_cols
    end % m = 1:num_rows
 
    errorNormPercentage = (sqrt(differenceSum) / sqrt(referenceSum)) * 100;