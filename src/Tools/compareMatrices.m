function compareMatrices(Const, matA, matB)
    %compareMatrices(Const)
    %   Usage:
    %       compareMatrices(Const, ZmatA, ZmatB)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containging settings of which solver to run
    %       ZmatA, ZmatB
    %           Two mxn matrices that will be compared
    %
    %   Output Arguments:
    %      None    
    %
    %   Description:
    %       Compares two m x n matrices with each other and displays
    %       the difference

    narginchk(3,3);
    message_fc(Const,sprintf('  Comparing matrices'));
    
    if (size(matA)~=size(matB))
        message(Const,sprintf('Cannot compare matrices of different sizes'));
        error('Cannot compare matrices of different sizes');
    end%if
        
    if(true)
        for mm = 1:1%size(matA,1)
            for nn = 1:1%size(matA,2)
                fprintf("matA(%d,%d) = %.5f + %.5f\n", mm,nn,real(matA(mm,nn)), imag(matA(mm,nn)));
                fprintf("matB(%d,%d) = %.5f + %.5f\n", mm,nn,real(matB(mm,nn)), imag(matB(mm,nn)));
            end
        end
    end
    
    %error('wag hier');
    
    difference = abs(matA-matB);
    figure    
    imagesc(difference)