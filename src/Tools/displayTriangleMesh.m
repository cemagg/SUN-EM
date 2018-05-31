function displayTriangleMesh(Const, triangle_corners, triangle_centres, nodes)
    %displayTriangleMesh
    %   Date: 2018.05.18
    %   Usage:
    %       displayTriangleMesh(Const, triangleData)
    %
    %   Input Arguments:
    %       Const: A global struct containing:
    %       triangle_corners:
    %           A (num_triangles x 3) array containing the nodal indices of the 3 corners
    %       triangle_centres:
    %           A (num_triangles x 3) array containing the X,Y and Z co-ordinate of the 
    %           centrepoint of each triangle
    %
    %   Output Arguments:
    %       None
    %
    %   Description:
    %       Plots a triangle mesh and also labels each triangle ID at the centre-point of each
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.18
    %   Stellenbosch University
    %   Email: dludick.sun.ac.za

    figure

    hold on;
    box on;
    grid on;

    number_of_triangles = length(triangle_corners(:,1));

    for i = 1:number_of_triangles
   
        % Extract the node indices
        v1 = triangle_corners(i,1);
        v2 = triangle_corners(i,2);
        v3 = triangle_corners(i,3);
    
        v1XYZ = [nodes(v1,1), nodes(v1,2), nodes(v1,3)];
        v2XYZ = [nodes(v2,1), nodes(v2,2), nodes(v2,3)];
        v3XYZ = [nodes(v3,1), nodes(v3,2), nodes(v3,3)];

        triangleXcoordinates = [v1XYZ(1) v2XYZ(1) v3XYZ(1) v1XYZ(1)];
        triangleYcoordinates = [v1XYZ(2) v2XYZ(2) v3XYZ(2) v1XYZ(2)];
        triangleZcoordinates = [v1XYZ(3) v2XYZ(3) v3XYZ(3) v1XYZ(3)];

        % Display the triangle
        line(triangleXcoordinates, triangleYcoordinates, triangleZcoordinates,...
            'Color', 'k', 'LineWidth', 2);

        % Display extra information when in debugging mode
        if (true)
    
            % Number the triangles in the centres
            tCentre = triangle_centres(i,:);
            text(tCentre(1),tCentre(2),tCentre(3),num2str(i),'FontSize',16);
        
            % Display the node numbers on the triangle
            text(v1XYZ(1),v1XYZ(2),v1XYZ(3),num2str(v1),'FontSize',16);
            text(v2XYZ(1),v2XYZ(2),v2XYZ(3),num2str(v2),'FontSize',16);
            text(v3XYZ(1),v3XYZ(2),v3XYZ(3),num2str(v3),'FontSize',16);
        end
    end

    % TO-DO: Display also the basis functions on the edges

    xlabel('X','FontSize',16);
    ylabel('Y','FontSize',16);
    zlabel('Z','FontSize',16);
