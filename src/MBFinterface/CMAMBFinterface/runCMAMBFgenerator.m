function [CMAMBF] = runCMAMBFgenerator(Const, Solver_setup, zMatrices, yVectors, xVectors)
    Z = zMatrices.values(1:Solver_setup.mom_basis_functions_per_array_element , 1:Solver_setup.mom_basis_functions_per_array_element);
    V = yVectors.values(1:Solver_setup.mom_basis_functions_per_array_element);
    iFek = runCMAMOMsolver(Z, V);%xVectors.Isol(1:Solver_setup.mom_basis_functions_per_array_element);
    R = -1.*real(Z);
    X = -1.*imag(Z);
    D = R\X;
    %[J, l] = eigs(D);
    %[J, l] =  eig(X,R,'chol');
    [J, l] =  eigs(D, 5, 'SM');
    s = length(l);
    l = diag(l);
    NumArrayEls = Solver_setup.num_finite_array_elements;
    CMAMBF.eigenvalues = l;
    CMAMBF.CMs = J;
    CMAMBF.alpha = zeros(length(l), 1);
    CMAMBF.sol = zeros(length(J), length(l)+1);
    
    iTot = 0;
    for n=1:size(l)
        Jn = J(:, n);
        ln = l(n);
        length(Jn);
        length(V);
        top = dot(Jn,V); 
        CMAMBF.alpha(n, 1) = top/(1+ln*1i);
        CMAMBF.sol(:, n) = CMAMBF.alpha(n, 1)*Jn;
        iTot = iTot + CMAMBF.sol(:, n);
    end
    sm = iFek - iTot;
    CMAMBF.sol(:, length(l)+1) = sm;

    redSolu = reduceCMAMBFs(Const, CMAMBF.sol);
    CMAMBF.redSol = redSolu;
    CMAMBF.numRedSol = length(redSolu(1,:));
end

