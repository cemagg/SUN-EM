% =========================================================================
% Nane: circular_array_driver.m
% Description: Creates a regular / random circular array
% =========================================================================
close all
clear all

clf
axis equal
hold on

% Some details of the freq., wavelength, etc.
c0 = 299792458;    % Speed of light in free-space
freq = c0;         % The operating frequency (70MHz (70e6) for SKA, c0 for dipoles)
lambda = c0/freq;  % The operating wavelength
k0 = 2*pi/lambda;  % The free-space wave-number
DEG2RAD = pi/180;  % For converting degrees to radians
RAD2DEG = 180/pi;  % For converting radians to degrees
randomise = true; % TO-DO: This is not yet working for the different number of elements / sector
vary_num_phi_elements_with_r = true;

% Define the number of elements in the radial direction (r)
num_r_elements = 3;
r_offset = 2; % The radial offset

% Define the number of elements in the phi direction
if (vary_num_phi_elements_with_r)
    num_phi_elements = zeros(num_r_elements,1);
    phi_off_deg = zeros(num_r_elements,1);
    
    num_phi_elements(1) = 6;
    num_phi_elements(2) = 12;
    num_phi_elements(3) = 24;
    
    for r_index = 1:num_r_elements
        phi_off_deg(r_index) = 360.0 / num_phi_elements(r_index);
    end
else
    num_phi_elements = 8;
    phi_off_deg = 360.0 / num_phi_elements;
end

% Total number of array elements
if (vary_num_phi_elements_with_r)
    total_elements = sum(num_phi_elements);
else    
    total_elements = num_r_elements * num_phi_elements;
end

% if (randomise)
%     % Calculate a random factor for the phi offsets (for each element) ..
%     % use only the first element's phi spacing
%     phi_off_deg_random_factors = phi_off_deg*0.1*rand(total_elements,1); % 10% deviation of element position
% end

% Loop now over all the elements and calculate the x and y co-ordinates
% (Possibly adding some random spacing in the phi direction)
element_xy = zeros(total_elements,2);
element_index = 0;
for r_index = 1:num_r_elements
    if (vary_num_phi_elements_with_r)
        num_phi_elements_total = num_phi_elements(r_index);
        phi_off_deg_local = phi_off_deg(r_index);
    else
        num_phi_elements_total = num_phi_elements;
        phi_off_deg_local = phi_off_deg;
    end
    for phi_index = 1:num_phi_elements_total
        % Increment the element index
        element_index = element_index + 1;                
        
        % Calculate the radial offset for this section
        if (~randomise)
            phi = (phi_index - 1)*phi_off_deg_local;
            r = r_offset + (r_index - 1)*r_offset;
        else
            % Calculate a random factor for the phi offsets (for each element) ..
            % use only the first element's phi spacing
            phi_off_deg_random_factor = phi_off_deg_local*0.03*rand(1,1);
            r_random_factor = r_offset*0.4*rand(1,1);
            
            phi = (phi_index - 1)*(phi_off_deg_local + phi_off_deg_random_factor);
            r = r_offset + (r_index - 1)*(r_offset + r_random_factor);
        end
        
        % Calculate the (x,y) position for the element
        element_xy(element_index,1) = r*cos(DEG2RAD*phi);
        element_xy(element_index,2) = r*sin(DEG2RAD*phi);
        
        plot(element_xy(element_index,1),element_xy(element_index,2),'ro', 'MarkerSize', 16);
        text(element_xy(element_index,1),element_xy(element_index,2),num2str(element_index),'FontSize',16);
        
    end % r_index = 1:num_r_elements
end%for

% Write now these finite array positions to a file
fid = fopen('array_layout.xml', 'w+');
% First write out the header information
fprintf(fid,'<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid,'<ArrayDistributionMatrix xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xsd="http://www.w3.org/2001/XMLSchema" Version="2.0.0.0" DesignFrequency="%d" FrequencyUnit="Hz" CoordinateUnit="Meters" xmlns="http://www.antennamagus.com/schemas/ArrayDistributionMatrix.xsd"> <Elements>\n',freq);

for element_index = 1:total_elements
    fprintf(fid,'            <Element Name="%d" X="%f" Y="%f" Z="0" Magnitude="1.0" Phase="%f" PhiRotation="0" ThetaRotation="0" GammaRotation="0" />\n',element_index,element_xy(element_index,1),element_xy(element_index,2),0.0);
end

fprintf(fid,'    </Elements>\n');
fprintf(fid,'</ArrayDistributionMatrix>');
    
% Close again the file
fclose(fid);