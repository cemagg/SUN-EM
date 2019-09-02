function [Zact] = runCMADGFMsolver(Const, Solver_setup, zMatrices, xVectors)
%RUNCMADGFMSOLVER Summary of this function goes here
%   Detailed explanation goes here
    numArrayEls = Solver_setup.num_finite_array_elements;
    Ndgfm = Solver_setup.mom_basis_functions_per_array_element;
    Zact = complex(zeros(Ndgfm,Ndgfm));
    Ndgfm = Solver_setup.mom_basis_functions_per_array_element;
    alphanm = complex(ones(Ndgfm, numArrayEls));
    for m=1:numArrayEls
        domain_basis_functions = Solver_setup.rwg_basis_functions_domains{m};
        alphanm(:,m) = (xVectors.Isol(domain_basis_functions,1));
    end
    for m=1:numArrayEls
        domain_m_basis_functions = Solver_setup.rwg_basis_functions_domains{m};

        [Zact, Uact, Vact] = calcZmn(Const,zMatrices,1,m,m,domain_m_basis_functions,domain_m_basis_functions);
        for n=1:numArrayEls
            if(m==n)
                continue;
            end
            domain_n_basis_functions = Solver_setup.rwg_basis_functions_domains{n};
            [Zmn, Umn, Vmn] = calcZmn(Const,zMatrices,1,m,n,domain_m_basis_functions, domain_n_basis_functions);
            Zact = Zact + Zmn;
            %Zact = Zact + alphanm(:, n).*Zmn;
        end
    end
end

