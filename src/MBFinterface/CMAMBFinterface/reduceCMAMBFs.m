function [reducedMBF] = reduceMBFJacques(Const, fullMBF)
%REDUCEMBFJACQUES Summary of this function goes here
%   Detailed explanation goes here
[U,S,V] = svd(fullMBF, 0);
thr = 10000000;

if size(S,1)>1
    numBases = sum(diag(S)./max(diag(S))>(1/thr),1);
    reducedMBF = U(:, 1:numBases);
else
    reducedMBF = U(:, 1:1);
end
end

