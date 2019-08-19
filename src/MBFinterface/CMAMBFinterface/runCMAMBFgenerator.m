function [CMAMBF] = runCMAMBFgenerator(Const, zMatrices, yVectors)
    Z = zMatrices.values;
    R = -1.*real(Z);
    X = -1.*imag(Z);
    D = R\X;
    %[J, l] = eigs(D);
    %[J, l] =  eig(X,R,'chol');
    [J, l] =  eigs(D, 27, 'SM');
    s = length(l);
    l = diag(l);
    CMAMBF.eigenvalues = l;
    CMAMBF.CMs = J;
    CMAMBF.alpha = zeros(length(l), 1);
    CMAMBF.sol = zeros(length(J), length(l));
    V = yVectors.values;
    for n=1:size(l)
        Jn = J(:, n);
        ln = l(n);
        length(Jn);
        length(V);
        top = dot(Jn,V); 
        CMAMBF.alpha(n, 1) = top/(1+ln*1i);
        CMAMBF.sol(:, n) = CMAMBF.alpha(n, 1)*Jn;
    end
    redSolu = reduceCMAMBFs(Const, CMAMBF.sol);
    CMAMBF.redSol = redSolu;
end

