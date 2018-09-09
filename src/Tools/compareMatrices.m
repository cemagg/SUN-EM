function compareMatrices(Const, refmat, mat)
    %compareMatrices(Const)
    %   Usage:
    %       compareMatrices(Const, refmat, mat)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containging settings of which solver to run
    %       refmat, mat
    %           Two mxn matrices that will be compared. The reference matrix is taken
    %           as refmat.
    %
    %   Output Arguments:
    %       None
    %
    %   Description:
    %       Compares two m x n matrices with each other by calculating the Matrix error 
    %       norm percentage. We also plot the difference between the matrix elements for
    %       visual comparison.

    narginchk(3,3);
    message_fc(Const,sprintf('Comparing matrices'));
    
    if (size(refmat)~=size(mat))
        message_fc(Const,sprintf('Cannot compare matrices of different sizes'));
        error('Cannot compare matrices of different sizes');
    end%if
       
    % Display a subset of the matrices 
    if(false)
        for mm = 1:100%size(matA,1)
            for nn = 1:100%size(matA,2)
                fprintf('refmat(%d,%d) = %.5f + %.5f\n', mm,nn,real(refmat(mm,nn)), imag(mat(mm,nn)));
                fprintf('mat(%d,%d)    = %.5f + %.5f\n', mm,nn,real(mat(mm,nn)), imag(mat(mm,nn)));
                fprintf('\n');
            end
        end
    end
    
    % Calculate the relative error norm %.
    Zerr = calculateMatrixErrorNormPercentage(refmat, mat);
    message_fc(Const,sprintf('  Rel. error norm. for Z compared to FEKO sol. %f percent',Zerr));    

    difference = abs(refmat-mat);        
    figure    
    imagesc(difference)

    % Also calcualte the maximum difference for the matrix elements
    matrix_diff = max(max(refmat-mat));
    message_fc(Const,sprintf('  Maximum difference in Z-matrix entry: %f + j*%f',real(matrix_diff),imag(matrix_diff)));    