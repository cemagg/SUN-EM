function displayTriangleMesh(Const, Solver_setup)
    %displayTriangleMesh
    %   Date: 2018.05.18
    %   Usage:
    %       displayTriangleMesh(Const, triangleData)
    %
    %   Input Arguments:
    %       Const: 
    %           A global struct with program flow settings
    %       Solver_setup: 
    %           The solver setup (basis geometry, basis functions, etc)
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
    
    % Extract the necessary data
    triangle_corners = Solver_setup.triangle_vertices;
    triangle_centres = Solver_setup.triangle_centre_point;
    nodes = Solver_setup.nodes_xyz;
    
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
    
    % Display the basis functions (also distinguish in which domains they
    % are). Note, currently this is only for domain decomposition. Should
    % be trivial to extend for general case. To show all basis functions, 
    % set plot_domain = 0;
    plot_domain = 0;
    % Only plot the internal basis functions associated with the domain
    internal_bfs_only = false;

    % We can also rather plot the generating sub-array RWGs
    plot_generating_sub_array_basis_functions = true;

    if (Solver_setup.disconnected_domains && plot_generating_sub_array_basis_functions)
        plot_generating_sub_array_basis_functions = false;        
    end %if

    if (~plot_generating_sub_array_basis_functions)
    
        % Plot for each of the domains the basis functions. We can also
        % display only the basis functions of a specific domain if 
        % domain_id <> -1 as set above
        for domain_index = 1:Solver_setup.number_of_domains
            if ((domain_index == plot_domain) || (plot_domain == 0))
                if (internal_bfs_only)                
                    domain_basis_functions = Solver_setup.rwg_basis_functions_internal_domains{domain_index};
                else
                    domain_basis_functions = Solver_setup.rwg_basis_functions_domains{domain_index};
                end
                
                for i = 1:length(domain_basis_functions)
                    bf_index = domain_basis_functions(i);
                    % Extract the centrepoint of the shared edge - this is where we
                    % will put the label
                    rwg_pC = Solver_setup.rwg_basis_functions_shared_edge_centre(bf_index,:);

                    text(rwg_pC(1),rwg_pC(2),rwg_pC(3),num2str(bf_index),...
                        'FontSize',16,'EdgeColor','red');
                    % TO-DO: Draw an arrow between the centre-point of triangle + and -
                end  %for
            end %if ((domain_index == plot_domain) || (plot_domain == 0))
        end %for domain_index = 1:Solver_setup.number_of_domains

    else
        % Plot for each of the generating sub-arrays the basis functions
        % (only applicable if we have connected domains). We can also
        % display only the basis functions of a specific domain if 
        % domain_id <> -1 as set above
        for ii = 1:Solver_setup.generating_subarrays.number_of_domains
            if ((ii == plot_domain) || (plot_domain == 0))
                
                % Extract the sub-array basis functions
                subarray_basis_functions = ...
                    Solver_setup.generating_subarrays.rwg_basis_functions_domains{ii};
                
                for jj = 1:length(subarray_basis_functions)
                    bf_index = subarray_basis_functions(jj);
                    % Extract the centrepoint of the shared edge - this is where we
                    % will put the label
                    rwg_pC = Solver_setup.rwg_basis_functions_shared_edge_centre(bf_index,:);

                    text(rwg_pC(1),rwg_pC(2),rwg_pC(3),num2str(bf_index),...
                        'FontSize',16,'EdgeColor','green');
                    % TO-DO: Draw an arrow between the centre-point of triangle + and -
                end  %for
            end %if ((domain_index == plot_domain) || (plot_domain == 0))
        end %for domain_index = 1:Solver_setup.number_of_domains
                        
    end%if

