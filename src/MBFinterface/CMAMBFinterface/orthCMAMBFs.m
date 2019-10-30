function [orthMBF] = orthCMAMBFs(fullMBF)
    %   Usage:
    %       redSolu = orthCMAMBFs(Const, CMAMBF.sol);
    %   Input arguments:
    %       Const
    %           A global struct, containing general data
    %       fullMBF
    %           Matrix with modal weighting coefficients as MBFs
    %   Output arguments:
    %       orthMBF
    %           Orthonormalized fullMBF matrix
    %   Description:
    %       This function orthonormalizes the generated MBFs
    %       the options I've been working on are the "orth()" function
    %       and the "svd()" function. When using the "svd()" function
    %       the U variable is used as the orth'd matrix. No reduction
    %       is being done, since CMA
   

    [U,S,V] = svd(fullMBF, 0);
    
    orthMBF = U; %orth(fullMBF);
    
end

