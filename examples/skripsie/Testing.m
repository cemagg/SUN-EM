
A = [1 2 3 1; 4 1 5 6; 7 1 8 9; 8 1 8 8; 9 10 11 12];

%A(:,:,2) = [10 11 12 2; 13 14 15 2; 1 6 17 18; 4 2 4 4];
%A(:,:,3) = [1 9 20 21; 2 2 23 24; 2 5 26 27; 4 1 5 6];
%A(:,:,4) = [1 9 20 21; 2 2 23 24; 2 5 26 27; 5 6 5 8];
row = 1;
for i = 1:2

    for j = 1:2
        for k = 1:4
            threeD(j,k,i) = A(row,k);

        end 
        row = row+1;
    end
end


