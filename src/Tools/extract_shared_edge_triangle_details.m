function [Solver_setup] = extract_shared_edge_triangle_details(Const, Solver_setup)
    %extract_shared_edge_triangle_details
    %   Date: 2018.05.23
    %   Usage:
    %       Solver_setup = extract_shared_edge_triangle_details(Const, Solver_setup)
    %
    %   Input Arguments:
    %       Const: 
    %           A global struct containing:
    %       Solver_setup
    %           A global struct containing solver setup details, e.g. frequency range,
    %           geometry data, basis function setup, etc.
    %   Output Arguments:
    %       Solver_setup:
    %           The above mentioned struct that is now updated with:
    %             Solver_setup.rwg_basis_functions_shared_edge_nodes
    %                 The two nodes defining the RWG shared edge vertex points
    %             Solver_setup.rwg_basis_functions_trianglePlusFreeVertex, Solver_setup.rwg_basis_functions_triangleMinusFreeVertex
    %                 The Tm+ and TM- free vertices
    %
    %   Description:
    %       Calculates some additional parameters associated with the RWG edge (e.g. shared edge nodal points, the Tm+ and Tm-
    %       free vertices).
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.253
    %   Stellenbosch University
    %   Email: dludick.sun.ac.za    

    narginchk(2,2);
    
    % Local debug flag to print some additional output
    LOCAL_DEBUG = false;

    % We also need to consider the number of domains - if present.
    % (typically if we are running the CBFM, DGFM, or other DD method)
    if (Const.domain_decomposition)
        % Domains are specified using labels in FEKO. The labels (as read
        % from the *.out file) also tells us how many domains are active
        Solver_setup.domain_labels = unique(Solver_setup.metallic_triangles_labels);
        Solver_setup.number_of_domains = size(Solver_setup.domain_labels,1);
        % Allocate now space to store the basis functions list of each
        % domain (will be used when extracting e.g. matrix entries)
        Solver_setup.rwg_basis_functions_domains = cell(Solver_setup.number_of_domains,1);
        
        % Initialise now each of the domain basis function lists
        for domain_index = 1:Solver_setup.number_of_domains
            Solver_setup.rwg_basis_functions_domains{domain_index} = [];
        end
        
        % Let's assume we have disjoint domains (flag will be updated below
        % if we find this is not the case)
        Solver_setup.disconnected_domains = true;
    end%if
    
    % Initialise the return values:
    Solver_setup.rwg_basis_functions_shared_edge_nodes = zeros(Solver_setup.num_metallic_edges,2);    
    Solver_setup.rwg_basis_functions_trianglePlusFreeVertex = zeros(Solver_setup.num_metallic_edges,1);
    Solver_setup.rwg_basis_functions_triangleMinusFreeVertex = zeros(Solver_setup.num_metallic_edges,1);
    Solver_setup.rwg_basis_functions_shared_edge_centre = zeros(Solver_setup.num_metallic_edges,3);

    % Loop over each of the shared edges and extract the additional details mentioned below:
    for edge_index = 1:Solver_setup.num_metallic_edges

        triangle_plus  = Solver_setup.rwg_basis_functions_trianglePlus(edge_index,1);
        triangle_minus = Solver_setup.rwg_basis_functions_triangleMinus(edge_index,1);
        
        triangle_plus_vertices = [Solver_setup.triangle_vertices(triangle_plus,1), ...
            Solver_setup.triangle_vertices(triangle_plus,2), ...
            Solver_setup.triangle_vertices(triangle_plus,3)];

        triangle_minus_vertices = [Solver_setup.triangle_vertices(triangle_minus,1), ...
            Solver_setup.triangle_vertices(triangle_minus,2), ...
            Solver_setup.triangle_vertices(triangle_minus,3)];

        % Find p1 and p2 (shared edge nodal endpoints - indices in the nodes array)
        endpoints = intersect(triangle_plus_vertices, triangle_minus_vertices);
        
        Solver_setup.rwg_basis_functions_shared_edge_nodes(edge_index,1) = endpoints(1);
        Solver_setup.rwg_basis_functions_shared_edge_nodes(edge_index,2) = endpoints(2);
        
        % Also calculate and store the centre-point of this node.
        p1XYZ = [Solver_setup.nodes_xyz(endpoints(1),1), ...
            Solver_setup.nodes_xyz(endpoints(1),2), ...
            Solver_setup.nodes_xyz(endpoints(1),3)];
        
        p2XYZ = [Solver_setup.nodes_xyz(endpoints(2),1), ...
            Solver_setup.nodes_xyz(endpoints(2),2), ...
            Solver_setup.nodes_xyz(endpoints(2),3)];
        
        %p2XYZ = [nodes(endpoints(2),1),nodes(endpoints(2),2),nodes(endpoints(2),3)];
        
        Solver_setup.rwg_basis_functions_shared_edge_centre(edge_index,:) = lineCentre(p1XYZ, p2XYZ);

        % Find free vertex associated with Tm+
        Solver_setup.rwg_basis_functions_trianglePlusFreeVertex(edge_index,1) = ...
            setdiff(triangle_plus_vertices, triangle_minus_vertices);

        % Find free vertex associated with Tm- (Note order of arguments in setdiff reversed now)
        Solver_setup.rwg_basis_functions_triangleMinusFreeVertex(edge_index,1) = ...
            setdiff(triangle_minus_vertices, triangle_plus_vertices);                

        % Associate this edge in a domain (or multiple domains if we have an interconnected problem)
        triangle_plus_label  = Solver_setup.metallic_triangles_labels{triangle_plus};
        triangle_minus_label = Solver_setup.metallic_triangles_labels{triangle_minus};
        
        % If we us domain decomposition, then we need to identify the basis
        % functions residing in each of the domains
        if (Const.domain_decomposition)
        
            % Locate now the domain index for triangle plus and minus
            domain_index_triangle_plus  = find(strcmp(Solver_setup.domain_labels,triangle_plus_label));
            domain_index_triangle_minus = find(strcmp(Solver_setup.domain_labels,triangle_minus_label));

            % Add the basis function now to this domain's list
            if (domain_index_triangle_plus == domain_index_triangle_minus)
                % Basis function entirely in domain. Add it to the list.
                Solver_setup.rwg_basis_functions_domains{domain_index_triangle_plus} = [ ...
                    Solver_setup.rwg_basis_functions_domains{domain_index_triangle_plus},
                    edge_index];
            else
                % Basis function on the boundary - we have detected
                % interconnected domains!
                Solver_setup.disconnected_domains = false;
                
                % Add the basis function to both domains.
                Solver_setup.rwg_basis_functions_domains{domain_index_triangle_plus} = [ ...
                    Solver_setup.rwg_basis_functions_domains{domain_index_triangle_plus},
                    edge_index];
                
                Solver_setup.rwg_basis_functions_domains{domain_index_triangle_minus} = [ ...
                    Solver_setup.rwg_basis_functions_domains{domain_index_triangle_minus},
                    edge_index];
            end%if

            if (LOCAL_DEBUG)
                fprintf(1,'Triangle %d (+) in domain %d\n',triangle_plus, domain_index_triangle_plus);
                fprintf(1,'Triangle %d (-) in domain %d\n',triangle_minus, domain_index_triangle_minus);
            end%if
            
        end%if (Const.domain_decomposition)
            
    end%for edge_index = 1:Solver_setup.num_metallic_edges

