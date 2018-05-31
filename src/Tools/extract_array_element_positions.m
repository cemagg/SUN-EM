function [x, y] = extract_array_element_positions(Const, xml_file)
    %extract_array_element_positions
    %   Date: 2018.03.24
    %   Usage:
    %       [x, y] = extract_array_element_positions(Const, 'array_pos.xml')
    %
    %   Input Arguments:
    %       Const
    %           A global struct containing program flow settings
    %       xml_file
    %           The *.xml file name that contains the positions of 
    %           the array elements (e.g. total_array_layout.xml)
    %   Output Arguments:
    %       x,y
    %           The x and y co-ordinates in [m]
    %
    %   Description: (Original, as taken from UCL)
    %      Extract the x and y co-ordinates from a XML file that is used to
    %      specify the finite array layout (as read by the 'FA' card in
    %      FEKO)
    %
    %   =======================
    %   Author: Danie Ludick
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za    

error(nargchk(2,2,nargin));
message_fc(Const,' ');

message_fc(Const,sprintf('    Reading array element positions form %s', xml_file));

% To read the XML file, we make use of a xml2struct script that is
% available on mathworks file exchange :
% https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct
struct = xml2struct(xml_file);

% Extract the number of array elements in the structure:
number_array_elements = length(struct.ArrayDistributionMatrix.Elements.Element);

% Allocate some space for the x and y co-ordinates:
x = zeros(1,number_array_elements);
y = zeros(1,number_array_elements);

% Loop over each of the element data in the struct and extract the x and y
% co-ordinates:
for index = 1:number_array_elements    
    x(1,index) = str2double(struct.ArrayDistributionMatrix.Elements.Element{index}.Attributes.X);
    y(1,index) = str2double(struct.ArrayDistributionMatrix.Elements.Element{index}.Attributes.Y);
    fprintf("      Extracting position of element %d = (%.2f, %.2f)\n",index,x(1,index),y(1,index));
end%for
