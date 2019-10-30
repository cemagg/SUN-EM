function [CMAMBF] = runCMA_MBFgenerator(Const, Solver_setup, zMatrices, yVectors, Modes)
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details    
    %       zMatrices
    %           The Z-matrices data
    %       yVectors
    %           The Yrhs-vector data
    %       yVectors
    %           The Xsol-vector data (i.e. MoM solution of FEKO)
    %       Modes
    %           Integer determining how many characteristic modes must be
    %           generated
    %   Output Arguments:
    %       CMAMBF
    %           Structs containing orthonormalized MBFs generated with
    %           Characteristing Modes, as well as how many MBFs has been
    %           used.
    %   Description:
    %       This function generates MBFs by using Characteristic Mode 
    %       Analysis, it also adds the option of adding a Source mode as an
    %       extra MBF, and the option of hybridizing the MBF generation
    %       with the DGFM
    %   
    %   The V variable is filled with the excitation values that aligns
    %   with the impedances of the first array element.
    
    numArrayEls = Solver_setup.num_finite_array_elements;
    numBasisFunc = Solver_setup.mom_basis_functions_per_array_element;
    %   IF statement for checking whether the DGFM hybridization should be
    %   done
    if(Const.sourceMode)
        CMAMBF.sol = complex(zeros(numArrayEls, numBasisFunc, Modes+1));
        CMAMBF.orthSol = complex(zeros(numArrayEls, numBasisFunc, Modes+1));
    else
        CMAMBF.sol = complex(zeros(numArrayEls, numBasisFunc, Modes));
        CMAMBF.orthSol = complex(zeros(numArrayEls, numBasisFunc, Modes));
    end
    
    CMAMBF.eigenvalues = zeros(numArrayEls, Modes);
    CMAMBF.CMs = zeros(numArrayEls, numBasisFunc, Modes);
    CMAMBF.alpha = zeros(numArrayEls, Modes);
    tic
    if(Const.DGFM) 
        %   The Z matrix is then filled with the DGFM function's return
        %   value which is Zact
        for i=1:numArrayEls
            Zi = zMatrices.values(Solver_setup.rwg_basis_functions_domains{i} , :);
            V = yVectors.values(Solver_setup.rwg_basis_functions_domains{i});
            Zact = runCMADGFMsolver(Solver_setup, Zi);
            R = -1.*real(Zact);
            X = -1.*imag(Zact);
            D = R\X;
            %   using the "eigs()" function to do eigen decomposition
            [J, eigVal] =  eigs(D, Modes, 'SM');
            %   Creating vector containing the eigenvalues from the matrix containing
            %   the eigenvalues on the diagonal of the matrix
            eigVal = diag(eigVal);
            CMAMBF.eigenvalues(i, :) = eigVal;
            CMAMBF.CMs(i, :, :) = J;
            %   This IF statement just checks whether to add a source mode as a MBF
            
            iTot = 0;
            %   This for loop calculates the modal weighting coefficients and
            %   multiplies it with the eigencurrents, the product of the two are
            %   then used as MBFs and added to the CMAMBF.sol
            for n=1:size(eigVal)
                Jn = J(:, n);
                eigValn = eigVal(n);
                numerator = dot(Jn,V); 
                CMAMBF.alpha(i,n) = numerator/(1+eigValn*1i);
                CMAMBF.sol(i, :, n) = CMAMBF.alpha(i, n)*Jn;
                iTot = iTot + CMAMBF.sol(i, :, n);
            end
            %   This IF statement checks whether the source mode should come from 
            %   subtracting iTot with the currents only on the element or on the
            %   DGFM Zact matrix element. It then adds the appropriate source mode
            if(Const.sourceMode)
                iFek = runCMAMOMsolver(Zact, V); 
                sm = iFek - transpose(iTot);
                CMAMBF.sm = sm;
                CMAMBF.sol(i, :, length(eigVal)+1) = sm;
            end
            %   Finally the MBFs are orthonormalized by the "orthCMAMBFs()"
            %   function
            
            orthSolu = orthCMAMBFs(squeeze(CMAMBF.sol(i, :, :)));
            CMAMBF.orthSol(i, :, :) = orthSolu;
            CMAMBF.numMBFs = length(orthSolu(1,:));
        end
    else
        for i=1:numArrayEls
            Z = zMatrices.values(Solver_setup.rwg_basis_functions_domains{i} , Solver_setup.rwg_basis_functions_domains{i});
            V = yVectors.values(Solver_setup.rwg_basis_functions_domains{i});
            
            R = -1.*real(Z);
            X = -1.*imag(Z);
            D = R\X;
            %   using the "eigs()" function to do eigen decomposition
            [J, eigVal] =  eigs(D, Modes, 'SM');
            %   Creating vector containing the eigenvalues from the matrix containing
            %   the eigenvalues on the diagonal of the matrix
            eigVal = diag(eigVal);
            CMAMBF.eigenvalues(i, :) = eigVal;
            CMAMBF.CMs(i, :, :) = J;
            %   This IF statement just checks whether to add a source mode as a MBF
            
            iTot = 0;
            %   This for loop calculates the modal weighting coefficients and
            %   multiplies it with the eigencurrents, the product of the two are
            %   then used as MBFs and added to the CMAMBF.sol
            for n=1:size(eigVal)
                Jn = J(:, n);
                eigValn = eigVal(n);
                numerator = dot(Jn,V); 
                CMAMBF.alpha(i,n) = numerator/(1+eigValn*1i);
                CMAMBF.sol(i, :, n) = CMAMBF.alpha(i, n)*Jn;
                iTot = iTot + CMAMBF.sol(i, :, n);
            end
            %   This IF statement checks whether the source mode should come from 
            %   subtracting iTot with the currents only on the element or on the
            %   DGFM Zact matrix element. It then adds the appropriate source mode
            if(Const.sourceMode)
                iFek = runCMAMOMsolver(Z, V); 
                sm = iFek - transpose(iTot);
                CMAMBF.sm = sm;
                CMAMBF.sol(i, :, length(eigVal)+1) = sm;
            end
            %   Finally the MBFs are orthonormalized by the "orthCMAMBFs()"
            %   function
            
            orthSolu = orthCMAMBFs(squeeze(CMAMBF.sol(i, :, :)));
            CMAMBF.orthSol(i, :, :) = orthSolu;
            CMAMBF.numMBFs = length(orthSolu(1,:));
        end
    end
    CMAMBF.generationTime = toc;
    %   Splitting the impedance matrix into real and imaginary parts and
    %   doing the necessary maths to prepare for eigen decomposition
    
end

