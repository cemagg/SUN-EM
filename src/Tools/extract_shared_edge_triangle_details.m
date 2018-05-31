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
    %             Solver_setup.rwg_basis_functions_sharededge_nodes
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

    error(nargchk(2,2,nargin));

    % Initialise the return values:
    Solver_setup.rwg_basis_functions_sharededge_nodes = zeros(Solver_setup.num_metallic_edges,2);    
    Solver_setup.rwg_basis_functions_trianglePlusFreeVertex = zeros(Solver_setup.num_metallic_edges,1);
    Solver_setup.rwg_basis_functions_triangleMinusFreeVertex = zeros(Solver_setup.num_metallic_edges,1);

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
        
        Solver_setup.rwg_basis_functions_sharededge_nodes(edge_index,1) = endpoints(1);
        Solver_setup.rwg_basis_functions_sharededge_nodes(edge_index,2) = endpoints(2);

        % Find free vertex associated with Tm+
        Solver_setup.rwg_basis_functions_trianglePlusFreeVertex(edge_index,1) = ...
            setdiff(triangle_plus_vertices, triangle_minus_vertices);

        % Find free vertex associated with Tm- (Note order of arguments in setdiff reversed now)
        Solver_setup.rwg_basis_functions_triangleMinusFreeVertex(edge_index,1) = ...
            setdiff(triangle_minus_vertices, triangle_plus_vertices);                

    end%for edge_index = 1:Solver_setup.num_metallic_edges

