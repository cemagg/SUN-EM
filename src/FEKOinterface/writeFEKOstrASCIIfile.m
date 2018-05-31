% Author: DJ Ludick
% Date: 02.23.2012
% Version: 1.0.0
% Description: Creates a *.str file containing the characteristic modes

% Parameters:
% Input: fileName     - the name of the ascii *.str file
%        numModes     - the Number of modes that are calculated
%                       ordered according to smallest magnitude of the
%                       corresponding eigenvalue
%        numBasisf   -  The number of RWG basis-functions of the problem
%        charModes   -  The (numBasisf x numBasisf) vector containing all 
%                       the eigen-currents. Note that only the first
%                       numModes are written to the out-file, i.e.
%                       'eigenModes.str'.

function create_CMmode_str_ascii_file(inFileName, numModes, numBasisf, charModes)

% Open the (Input) file
fidIn = fopen(inFileName, 'r');

% OutFile
outFileName = 'eigenModes.str';
fidOut = fopen(outFileName, 'w');

% Read the header information from the input file and write this to the
% output file
tline = fgets(fidIn);
fprintf(fidOut,tline);

tline = fgets(fidIn);
fprintf(fidOut,tline);

tline = fgets(fidIn);
num_file_coeff=sscanf(tline,'%d');
if (num_file_coeff ~= numBasisf)
    error ('File Read Error');
end 
fprintf(fidOut,tline);

tline = fgets(fidIn);
fprintf(fidOut,tline);

tline = fgets(fidIn);
fprintf(fidOut,tline);

% Now write the characteristic modes to the eigenModes.str file.
for ii=1:numModes,
    % Write out for each characteristic mode its value to the *.out
    % file
    for jj=1:numBasisf
        fprintf ( fidOut, '(%20.10e,%20.10e)\n', real(charModes(jj,ii)), imag(charModes(jj,ii)) );
    end    
% write a delimiter line that seperates the individual modes    
    fprintf(fidOut,'--- separation of different blocks ---\n');
end

% Close the file
fclose(fidOut);
fclose(fidIn);