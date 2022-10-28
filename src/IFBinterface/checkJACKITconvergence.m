function [converged] = checkJACKITconvergence(Const, zMatricesFEKO)
    %runJACKITsolver v0.1
    %   Date: 02.08.2013
    %   Usage:
    %       [converged] = checkJACKITconvergence(Const, zMatrices)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       zMatrices
    %           The Z-matrices data
    %
    %   Output Arguments:
    %       converged
    %           Logical variable indicating that the Jacobi Iterative
    %           solution will converge to the full-wave MoM solution
    %           (true) - or not (false)
    %
    %   Description:
    %       Checks whether the Jacobi Iterative Solver, as described in [1] and
    %       [2] and adapted for the DGFM in [3], will converge. According
    %       to Dr. Rob Maaskant in [3], Appendix A, pp. 9-10, convergence will be achieved
    %       if "the magnitude of the eigenvalue of the principal eigenvector of 
    %       (Zon)^(-1(Zoff) is less than unity."
    %
    %   Assumptions:
    %        - All domains are the same size (i.e. contains the same number of unknowns)
    %        - The principle eigenvalue is also the largest eigenvalue
    %
    %   References:
    %   [1] Y. Brand, A. K. Skrivervik, J. R. Mosig, F. E. Gardoil, "New
    %       iterative integral equation technique for multilayered printed
    %       array antennas," in Mathematical Methods in Electromagnetic Theory,
    %       Kharkov, Ukraine, Jun, 1998, pp. 615-617.
    %   [2] A. C. Polycarpou, "Evaluation of stationary block iterative
    %       techniques for the solution of finite arrays using fe-bi method
    %       and domain decomposition," in Proc. European Conference on
    %       Antennas and Propagation (EuCAP), Nice, France, Nov. 2006, 
    %       pp. 1-6
    %   [3] D.J. Ludick, R. Maaskant, et. al, "Efficient Analysis of Large 
    %       Irregular Antenna Arrays using the Domain Green's Function Method", 
    %       Special Issue on Antennas Propagation Submission, 2013 (pending)
    %   =======================
    %   Written by Danie Ludick on August 02, 2013.
    %   Last updated on June 24, 2013.
    %   EMSS-SA (Pty) Ltd
    %   Email: dludick@emss.co.za

    error(nargchk(2,2,nargin));
        
    if (~ Const.JACKITcheckConvergence) 
        % Return immediately if convergence check not required.
        converged = true;
        return;
    end%if
    
    message(Const,sprintf('  Checking convergence for Jacobi Iterative solver'));
    converged = true;
    
    % Initialise the return values
    Nmom = Const.numMoMbasis;
    Ndom = Const.numMoMbasisPerElement;
    numArrayEls = Const.numArrayElements;

    % Switch additional debug information for the Eigs routine on only in
    % DEBUG mode
    if  (Const.debug)
        OPTS.disp = 1;
    else
        OPTS.disp = 0;
    end
    
    % TO-DO: - This currently works for a single Yrhs and freq. -> repeat also for multiple Yrhs
    %        - Also assumed here are 1-to-1 mapping, i.e. Const.arrayMappingVector
    %          not yet used.
    Zon = complex(zeros(Ndom,Ndom));
    Zoff = complex(zeros(Ndom,Ndom));
    for m=1:numArrayEls        
        % Extract first the self-interaction / diagonal matrix here (Zon)
        Zon = zMatricesFEKO.values((m-1)*Ndom+1:m*Ndom,(m-1)*Ndom+1:m*Ndom);
        for n=1:numArrayEls
            % Extract the coupling matrix, Zoff
            Zoff = zMatricesFEKO.values((m-1)*Ndom+1:m*Ndom,(n-1)*Ndom+1:n*Ndom);

            % Calculate the principle eigenvalue of (Zon)^(-1) x Zoff
            OPTS.display = 0;
            [EigVec, EigVal, flag] = eigs(inv(Zon)*Zoff,[],1,'LM',OPTS);
            if (flag ~= 0)
                % Critical error - no convergence for Eigs
                error(sprintf('No convergence for Eigs() for Domain M=%d,N=%d',m,n));
                message(Const,sprintf('No convergence for Eigs() for Domain M=%d,N=%d',m,n));
            end%if
            if (Const.debug) 
                % Useful debug output
                message(Const,sprintf('  Principle |Eigenvalue| for (Zon(%d,%d))^(-1) x Zoff(%d,%d) = %f',m,m, m,n, abs(EigVal(1)) ));
            end%if
            if ((abs(EigVal(1)) - 1.0) > Const.eps)
                converged = false;
                %error(sprintf('No convergence for Domain M=%d,N=%d with EigVal: %f',m,n, abs(EigVal(1))));
                message(Const,sprintf('No convergence for Domain M=%d,N=%d with EigVal: %f',m,n, abs(EigVal(1))));
            end%if
        end
    end
        
    % Formating based on where the routine is called from
    message(Const,sprintf('  Finished Jacobi Iterative convergence check with converged: %d ',converged));