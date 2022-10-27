% Example to show how to calculate the frobenius norm associated with
% 2 x complex matrices, where 1 is the reference solution and the other
% one that we want to figure out what the accuracy is

% Generate some reference data.
refMatrix = complex(zeros(3,3));
refMatrix = refMatrix + rand(3) + 1i.*rand(3);

% First step is to check that we get 0% error when we pass our refData to the function
errorNormPercentage = calculateMatrixErrorNormPercentage(refMatrix, refMatrix);
fprintf(1,"Frobenius Error Norm Percentage = %f percent\n",errorNormPercentage);

% Now generate another matrix (e.g. with another technique used to
% calculate the matrix elements). Here we just use the same values as in
% refMatrix, but we add some small random numbers (you can adjust the
% accuracy by adjusting the division factor -> e.g. 1000 will result in a
% very small Frob norm (i.e. the matrix elements are essentailly the same -
% verify this by looking at the elements).
newMatrix = refMatrix + rand(3)./10 + 1i.*rand(3)./10;

% Calculate the error norm percentage of the new matrix by comparing it
% against the reference data.
errorNormPercentage = calculateMatrixErrorNormPercentage(refMatrix, newMatrix);
fprintf(1,"Frobenius Error Norm Percentage = %f percent\n",errorNormPercentage);
