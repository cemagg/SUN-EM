% TO-DO: This script needs to be converted to a function, that
%        takes parameters and returns the file-pointer / location
%        to the array_layout.xml file that is read by FEKO.
%        Also the comments need to be updated, in order to be consistent
%        with that of FEKO.

%  Generates a set of random points, with a minimum average distance
%  between them
clear 

clf
axis equal
hold on

% Some details of the freq., wavelength, etc.
c0 = 299792458;    % Speed of light in free-space
freq = c0;       % The operating frequency (70MHz (70e6) for SKA, c0 for dipoles)
lambda = c0/freq   % The operating wavelength
k0 = 2*pi/lambda;  % The free-space wave-number
DEG2RAD = pi/180;  % For converting degrees to radians
RAD2DEG = 180/pi;  % For converting radians to degrees

START_SCHEME = 1; % How do we start with the list of elements that we will include
                  % 0 - old technique, start from element with most
                  %     potential candidates
                  % 1 - new technique, start from element closest to centre
                  
FA_OUTPUT = 1     % How we write the array layout to the file
                  % 0 - standard 'FA' card format - distribution follows in
                  %     *.pre file (has to be copied manually from the
                  %     file and pasted in the *.pre file)
                  % 1 - write out a *.xml MAGUS finite array distribution
                  % 2 - write out a format used in many FEKDDM-x projects:
                  %     #arr_posX[1] = 0.5*#lam        ** Element 1 X position
                  %     #arr_posY[1] = 0.75*#lam       ** Element 1 Y position
                  %     #arr_posZ[1] = 0               ** Element 1 Z position

% See issue FEKDDM-5.2: Updated now this script to write out only an
% array_layout.xml file, corresponding to the first quadrant (we will
% use this in FEKO to create a circular array for the DGFM)
FIRST_QUADRANT_ONLY = 0;
                  
x1=0;                  % The unit circle X co-ordinate
y1=0;                  % The unit circle X co-ordinate
rc=25*lambda;          % The unit circle size (for large arrays, this should be ~25xlambda +)

numTotPoints = 100     % The total number of points that will be calculated
distAvMin = lambda*2 % The minimum average distance between the samples (make 1.2 for SKA)

% declare now some variables here for the X, Y co-ordinates and also the
% distance to each (R)
X = zeros(numTotPoints,1);
Y = zeros(numTotPoints,1);
R = zeros(numTotPoints,1);

[x,y,z] = cylinder([rc rc],200);
plot(x(1,:)+x1,y(1,:)+y1,'r')

for t=1:numTotPoints %loop until doing numTotPoints points inside the circle
[x, y]=cirrdnPJ(x1,y1,rc);
plot(x,y,'x')
% Display the node numbers on the triangle
%text(x,y,num2str(t),'FontSize',16);
%pause(0.01) %if you want to see the point being drawn
% Store now the X and Y co-cordinate
X(t) = x;
Y(t) = y;
R(t) = sqrt(x^2 + y^2);
end

% Loop now over the points and calculate the distance between each and
% store the result as a 2-dimensional array containing therefore 0s on
% diagonal
distArray = zeros(numTotPoints, numTotPoints);
for m=1:numTotPoints
    for n=1:numTotPoints
        distArray(m,n) = sqrt( (X(m)-X(n))^2 + (Y(m)-Y(n))^2 );
    end
end
% Print out the distance array
distArray;

% Now loop over all the elements in the distance array and keep all elements
% for which the distance from that element to all the others are more or equal 
% to the distAvMin (note, this is only for the case where m <> n). The IDs of
% the points that fulfills this requirement is stored in elemIDs.

elemIDs = zeros(numTotPoints, numTotPoints);
%elemIDs(:) = -1; % Initialise all the elemIDs to -1
for m=1:numTotPoints
    for n=m+1:numTotPoints
        if (distArray(m,n) >= distAvMin) 
            elemIDs(m,n) = 1;
            elemIDs(n,m) = 1;
        end
        % set the diagonal element to be 1        
    end
    elemIDs(m,m) = 1;
end

% Print out the element IDs
elemIDs;

% Find an element to start from
if (START_SCHEME == 0)
    % Loop over all the rows in the matrix and find the one with the most
    % non-zero entries (this will be a good starting point combination for the
    % candidate elements)
    numNonZeroEntries = 0;
    for m=1:numTotPoints
        if (length(find(elemIDs(m,:))) >= numNonZeroEntries)
            candidateRow = m;
            numNonZeroEntries = length(find(elemIDs(m,:)));
            candidateEntries = find(elemIDs(m,:));
        end
    end
else
    % Start from the element that is closest to the centre of the circle
    [minR,candidateRow] = min(R)
    numNonZeroEntries = length(find(elemIDs(m,:)));
    candidateEntries = find(elemIDs(m,:));
end

% Write out the data obtained above, i.e. the candidateRow index, the
% number of non-zero entries, as well as the element positions
candidateRow;
numNonZeroEntries;
candidateEntries;

% Now that we have the candidateRowEntries, we need to check the distances
% between all of them, and discard those that are not sufficient.
for m=1:numNonZeroEntries
    % Check whether we have already discarded this entry from the list
    % if so, then we just continue to the next one
    if (candidateEntries(m)<=0)
        continue;
    end 
    for n=1:numNonZeroEntries        
        % Check whether we have already discarded this entry from the list
        % if so, then we just continue to the next one
        if (candidateEntries(n)<=0)
            continue;
        end 
        if (elemIDs(candidateEntries(n),candidateEntries(m)) ~= 1)
        % Cannot use this entry - mark as -1 and continue with next one
            candidateEntries(m) = 0;
            break;
        end
    end
end

candidateEntries;
% Extract now only the final entries
finalEntries = candidateEntries(find(candidateEntries > 0));
numFinalEntries = length(finalEntries)

% 28.11.2012:
% ----------
% We now have the final entries, but we need to add the correct phasing.
% Currently, we only allow a phasing along the Z/X-axis for a planar array
% located on the X-Y plane:
elPhase = zeros(numFinalEntries,1);

%before we continue, the elements need to be sorted in ascending X-values
[Xsorted perm] = sort(X(finalEntries(:)));
% Also permute now the finalEntries to be in this order
Ysorted = Y(finalEntries(perm));

% Plot the final points
for el = 1:numFinalEntries
    plot(Xsorted(el),Ysorted(el), 'ro', 'MarkerSize', 16)
    text(Xsorted(el),Ysorted(el),num2str(el),'FontSize',16);
    xyString = sprintf('(%5.2f,%5.2f)',Xsorted(el),Ysorted(el));
    text(Xsorted(el)+0.25,Ysorted(el),xyString,'FontSize',10);
end

% Now calculate for each of the elements along the X-axis, the correct
% phase - use a modified version of Eq. (8.7) in "Computational Electromagnetics 
% for RF and Microwave Engineering", DB Davidson. 

% Scan angle, Theta, in radians
Theta = 0 * DEG2RAD
% Calculate now the correct phasing of each entry
phase = 0
for el = 1:numFinalEntries
    % First calculate however the progressive phase constant
    % Calculate the distance, relative to the last element
    if (el == 1)
       dist = 0; 
    else
       dist = abs(Xsorted(el) - Xsorted(el-1))
    end
    delfz = -1 * k0 * dist * sin(Theta); % Multiply here the phase offset by -1, 
                                         % for thescan angle to be in +Theta direction
    % Now phase the elements correctly
    phase = phase + delfz * RAD2DEG;
    elPhase(el) = phase; 
end

% Now that we have the final entries - write out the X and Y co-ordinates
% in  a format that can be used with the 'FA'-card (either directly in 'FA'
% card output format that can be read in PREFEKO, or as an XML file that
% can also be read by the 'FA' card
if (FA_OUTPUT == 0)
    fid = fopen('arrayXYcoordinates.txt', 'w+');
    % Write the 'FA' card (complete with the number of samples that were
    % calculated.
    % TO-DO: Danie, element phase has not already been added here
    fprintf(fid,'FA: 1 : 2 : %d :  :  :  :  :  :  :  : 0\n',numFinalEntries);
    for el = 1:numFinalEntries    
        fprintf(fid,'  :  :  :  :  :  : %f : %f : 0 :  :  : 0 : 0 : 0\n',Xsorted(el),Ysorted(el));
    end
elseif (FA_OUTPUT == 2)
    fid = fopen('arrayXYcoordinates.txt', 'w+');
    % Consistent with some of the FEKDDM-x projects for the *.pre files -
    % see e.g. FEKDDM-5.1.4: EEPs for strip dipole arrays
    fprintf(fid,'#arr_numEls = %d ** The number of array elements\n',numFinalEntries);
    for el = 1:numFinalEntries    
        fprintf(fid,'#arr_posX[%d] = %f   ** Element %d X position\n',el,Xsorted(el),el);
        fprintf(fid,'#arr_posY[%d] = %f   ** Element %d Y position\n',el,Ysorted(el),el);
        fprintf(fid,'#arr_posZ[%d] = %f   ** Element %d Z position\n',el,0.0,el);
        fprintf(fid,'#arr_rotX[%d] = %f   ** Element %d X rotation\n',el,0.0,el);
        fprintf(fid,'#arr_rotY[%d] = %f   ** Element %d Y rotation\n',el,0.0,el);
        fprintf(fid,'#arr_rotZ[%d] = %f   ** Element %d Z rotation\n',el,0.0,el);
        fprintf(fid,'\n');
    end
    
else % Standard array_layout.xml format    
    fid = fopen('array_layout.xml', 'w+');
    % First write out the header information
    fprintf(fid,'<?xml version="1.0" encoding="utf-8"?>\n');
    fprintf(fid,'<ArrayDistributionMatrix xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xsd="http://www.w3.org/2001/XMLSchema" Version="2.0.0.0" DesignFrequency="%d" FrequencyUnit="Hz" CoordinateUnit="Meters" xmlns="http://www.antennamagus.com/schemas/ArrayDistributionMatrix.xsd"> <Elements>\n',freq);
    % See issue FEKDDM-5.2: Write out only the first quadrant of elements if FIRST_QUADRANT_ONLY is set to 1
    if (FIRST_QUADRANT_ONLY == 0)
        % Full 360 deg. circle
        for el = 1:numFinalEntries
            fprintf(fid,'            <Element Name="%d" X="%f" Y="%f" Z="0" Magnitude="1.0" Phase="%f" PhiRotation="0" ThetaRotation="0" GammaRotation="0" />\n',el,Xsorted(el),Ysorted(el),elPhase(el));
        end
    else
        % First quadrant, i.e. 0 deg. <= theta <= 90 deg.
        for el = 1:numFinalEntries
            if ( (Xsorted(el) > 0) && (Ysorted(el) > 0) )
                fprintf(fid,'            <Element Name="%d" X="%f" Y="%f" Z="0" Magnitude="1.0" Phase="%f" PhiRotation="0" ThetaRotation="0" GammaRotation="0" />\n',el,Xsorted(el),Ysorted(el),elPhase(el));
            end                
        end
    end
    fprintf(fid,'    </Elements>\n');
    fprintf(fid,'</ArrayDistributionMatrix>');
end

% Close again the file
fclose(fid);