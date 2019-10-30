function [Zact] = runCMADGFMsolver(Solver_setup, Zi)
%RUNCMADGFMSOLVER Summary of this function goes here
%   Detailed explanation goes here
    numArrayEls = Solver_setup.num_finite_array_elements;
    Ndgfm = Solver_setup.mom_basis_functions_per_array_element;
    Zact = complex(zeros(Ndgfm,Ndgfm));
    alpha = 1;
    
    
    
    for n=1:numArrayEls
        domain_n_basis_functions = Solver_setup.rwg_basis_functions_domains{n};
        Zmn = Zi(:, domain_n_basis_functions);
        Zact = Zact + alpha*Zmn;
    end
end

