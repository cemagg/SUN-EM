% Scan angle, Theta, in radians
Theta = 60 * DEG2RAD
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
else
    fid = fopen('array_layout.xml', 'w+');
    % First write out the header information
    fprintf(fid,'<?xml version="1.0" encoding="utf-8"?>\n');
    fprintf(fid,'<ArrayDistributionMatrix xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xsd="http://www.w3.org/2001/XMLSchema" Version="2.0.0.0" DesignFrequency="%d" FrequencyUnit="Hz" CoordinateUnit="Meters" xmlns="http://www.antennamagus.com/schemas/ArrayDistributionMatrix.xsd"> <Elements>\n',freq);
    for el = 1:numFinalEntries
        % TO-DO: Danie, later we can add here the effect of phasing the
        % array
        fprintf(fid,'            <Element Name="%d" X="%f" Y="%f" Z="0" Magnitude="1.0" Phase="%f" PhiRotation="0" ThetaRotation="0" GammaRotation="0" />\n',el,Xsorted(el),Ysorted(el),elPhase(el));
    end
    fprintf(fid,'    </Elements>\n');
    fprintf(fid,'</ArrayDistributionMatrix>');
end

% Close again the file
fclose(fid);