function [macro] = runJacquesMBF(Const, Solver_setup, zMatrices, yVectors, xVectors)
%RUNJACQUESMBF Summary of this function goes here
%   Detailed explanation goes here
macro = [];
Nmom = Solver_setup.num_mom_basis_functions;
numSols = xVectors.numSols;
numArrayEls = Solver_setup.num_finite_array_elements;
max_primaries_per_domain = 1;
macro.RedIsol = complex(zeros(Nmom,numArrayEls,numArrayEls,numSols));
macro.PrimIsol = complex(zeros(Nmom, 1, numArrayEls, numSols));
macro.SecIsol = complex(zeros(Nmom,max_primaries_per_domain*(numArrayEls-1),numArrayEls,numSols));
macro.numPrimMBFs = zeros(numArrayEls,numSols); % Number of Prim. MBFs / solution config.
macro.numSecMBFs = zeros(numArrayEls,numSols); % Number of Sec.  MBFs / solution config.
domain_indices = Solver_setup.rwg_basis_functions_domains{1}; 
[L,U] = lu(zMatrices.values(domain_indices, domain_indices));

for solNum = 1:numSols
    macro.numPrimMBFs(:, solNum) = 0;
    for m=1:numArrayEls
        if (Const.is_array_element_active(m,solNum))
        
           domain_basis_functions = Solver_setup.rwg_basis_functions_domains{m};
           macro.numPrimMBFs(m,solNum) = macro.numPrimMBFs(m,solNum) + 1;
                    
           b = L\yVectors.values(domain_basis_functions,solNum);
           macro.PrimIsol(domain_basis_functions,1,m,solNum) = U\b;
        end
    end
    
    for m=1:numArrayEls
         domain_m_basis_functions = Solver_setup.rwg_basis_functions_domains{m};
         count = 0;
         for n=1:numArrayEls
             if (m == n)
                        %ignore self-coupling - primary MBF already calculated
                 continue;
             end
             if (Const.is_array_element_active(n,solNum))
                 count = count + 1;
                 macro.numSecMBFs(m,solNum) = count;  
                 domain_n_basis_functions = Solver_setup.rwg_basis_functions_domains{n};
                 Zcoupl = zMatrices.values(domain_m_basis_functions,domain_n_basis_functions);
                 Vcoupl =  - Zcoupl * macro.PrimIsol(domain_n_basis_functions,1,n,solNum);
                 b = L\Vcoupl;
                 macro.SecIsol(domain_m_basis_functions,count,m,solNum) = U\b;
                 
             end
         end
    end
    macro.totRedMBFs  = 0;
    if(Const.useMBFreduction)
        for m=1:numArrayEls
            MBF = [macro.PrimIsol(:,1:macro.numPrimMBFs(m,solNum),m,solNum) macro.SecIsol(:,1:macro.numSecMBFs(m,solNum),m,solNum)];
            redMBF = reduceMBFJacques(Const, MBF);
            macro.numRedMBFs(m,1) = size(redMBF,2);
            macro.RedIsol(:,1:size(redMBF,2),m,solNum) = redMBF;
            
        end
    end
    for n=1:numArrayEls
       macro.totRedMBFs = macro.totRedMBFs + macro.numRedMBFs(n,solNum);   
    end
end

