
ldamn = 5;

amn = zeros(ldamn,ldamn);

% Fill a 5x5 matrix
for m = 1:5
    for n = 1:5
        % Use 1-dimensional indexing to store the entry
        amn(m,n) = m + (n-1)*ldamn;
    end
end

% Output the matrix
amn

bfs = [1 2 4];
amn(bfs:bfs)
