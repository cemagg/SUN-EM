function [CMAMBF] = runCMA_MBFgenerator(Const, Solver_setup, zMatrices, yVectors, DGFM, SM, Modes)
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details    
    %       zMatrices
    %           The Z-matrices data
    %       yVectors
    %           The Yrhs-vector data
    %       xVectors
    %           The Xsol-vector data (i.e. MoM solution of FEKO)
    %       DGFM
    %           Boolean determining whether the DGFM hybridization will be
    %           used.
    %       SM
    %           Integer determining which source mode system will be used,
    %           or whether a source mode system will be used at all.
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
    V = yVectors.values(1:Solver_setup.mom_basis_functions_per_array_element);
    %   IF statement for checking whether the DGFM hybridization should be
    %   done
    if(DGFM) 
        %   The Z matrix is then filled with the DGFM function's return
        %   value which is Zact
        Z = runCMADGFMsolver(Const, Solver_setup, zMatrices);
    else
        %   The Z matrix is filled with the impedance values of the first
        %   array element
        Z = zMatrices.values(1:Solver_setup.mom_basis_functions_per_array_element , 1:Solver_setup.mom_basis_functions_per_array_element);
    end
    %   Splitting the impedance matrix into real and imaginary parts and
    %   doing the necessary maths to prepare for eigen decomposition
    R = -1.*real(Z);
    X = -1.*imag(Z);
    D = R\X;
    %   using the "eigs()" function to do eigen decomposition
    [J, l] =  eigs(D, Modes, 'SM');
    %   Creating vector containing the eigenvalues from the matrix containing
    %   the eigenvalues on the diagonal of the matrix
    l = diag(l);
    CMAMBF.eigenvalues = l;
    CMAMBF.CMs = J;
    CMAMBF.alpha = zeros(Modes, 1);
    %   This IF statement just checks whether to add a source mode as a MBF
    if(SM)
        CMAMBF.sol = zeros(length(J), Modes+1);
    else
        CMAMBF.sol = zeros(length(J), Modes);
    end
    iTot = 0;
    %   This for loop calculates the modal weighting coefficients and
    %   multiplies it with the eigencurrents, the product of the two are
    %   then used as MBFs and added to the CMAMBF.sol
    for n=1:size(l)
        Jn = J(:, n);
        ln = l(n);
        top = dot(Jn,V); 
        CMAMBF.alpha(n, 1) = top/(1+ln*1i);
        CMAMBF.sol(:, n) = CMAMBF.alpha(n, 1)*Jn;
        iTot = iTot + CMAMBF.sol(:, n);
    end
    %   This IF statement checks whether the source mode should come from 
    %   subtracting iTot with the currents only on the element or on the
    %   DGFM Zact matrix element. It then adds the appropriate source mode
    if(SM)
        iFek = runCMAMOMsolver(Z, V); 
        sm = iFek - iTot;
        CMAMBF.sm = sm;
        CMAMBF.sol(:, length(l)+1) = sm;
    elseif(SM == 2)
        Z = zMatrices.values(1:Solver_setup.mom_basis_functions_per_array_element , 1:Solver_setup.mom_basis_functions_per_array_element);
        iFek = runCMAMOMsolver(Z, V);
        sm = iFek - iTot;
        CMAMBF.sm = sm;
        CMAMBF.sol(:, length(l)+1) = sm;
    end
    %   Finally the MBFs are orthonormalized by the "orthCMAMBFs()"
    %   function
    orthSolu = orthCMAMBFs(Const, CMAMBF.sol);
    CMAMBF.orthSol = orthSolu;
    CMAMBF.numMBFs = length(orthSolu(1,:));
end

