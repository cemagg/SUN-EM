   function [I] = runCVXsolver(Const, Z, V)
    %runCVXsolver
    %   Usage:
    %       [I] = runCVXsolver(Const, Z, V)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Z
    %           The Z-matrix (typically and underdetermined system)
    %       V
    %           The RHS-vector data
    %
    %   Output Arguments:
    %       I
    %           argmin(||I||lo) subject to ZI=V, where ||I||lo measures the
    %           the number of non-zero entreis in I (i.e. the non-convex 
    %           optimisation problem.). See [1] for more detail. We essentially
    %           aim to solve the sparsest solution for I.
    %
    %   Description:
    %       Runs the CVX Matlab Software for Disciplined Convex Programming to 
    %       solve an undertermined system [2], following the method discussed in [1].
    %
    %   References:
    %   [1] R. Maaskant, D. Ludick, C. Bencivenni, M. V. Ivashina, D. Davidson,
    %       "A Compressed Sensing Technique Applied to the Characteristic Basis Function Method"
    %       2015 [to be published]
    %   [2] M. Grant and S. Boyd, "CVX Matlab software for disciplined convex programming, 
    %       version 2.1", http://www.cvxr.com/cvx, Mar. 2014

    error(nargchk(3,3,nargin));

    message(Const,sprintf('  Running CVX solver'));
    
    I_new=ones([size(Z,2) 1]); I_old=I_new*0;
    eps=1E-10;

    while norm(I_old-I_new)/norm(I_new)*100 > eps
        I_old=I_new;
        A=diag((abs(I_old)+eps).^(-1));
        
        n=size(I_new);
        
        cvx_clear;
        cvx_begin
        cvx_quiet(true);
        
        variable I(n) complex
        minimize( norm(A*I, 1 ) )
        subject to
        Z * I == V
        cvx_end
        
        if(strcmp(cvx_status,'Infeasible'))
            message(Const, 'ERROR: Infeasable solution');
            error('ERROR: Infeasable solution');
        elseif(strcmp(cvx_status,'Failed'))
            message(Const, 'ERROR: Failed');
            error('ERROR: Failed');
        end
        I_new=I;
    end