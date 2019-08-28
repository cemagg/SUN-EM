function [Imom] = runCMAMOMsolver(zMatrix,yVector)
    % LU-decomposition of the Z-matrix
    [L,U] = lu(zMatrix);
    b = L\yVector;
    Imom = U\b;
end

