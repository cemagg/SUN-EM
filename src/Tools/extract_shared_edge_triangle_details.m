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

        % Number of domains should correspond to the number of finite array
        % elements. TO-DO: Needs some refactoring when we allow different 
        % type of domains.
        Solver_setup.number_of_domains = Solver_setup.num_finite_array_elements;

        % Domains are specified using labels in FEKO. The labels (as read
        % from the *.out file) also tells us how many domains are active. First
        % read ALL the domain labels. Thereafter, we group all the labels ending
        % with @1, @2, etc. i.e. the unique way in which labels are prepended when 
        % geometry is arrayd using the 'FA' card.
        Solver_setup.all_domain_labels = unique(Solver_setup.metallic_triangles_labels);
        
        % There has to be equal number of labels / domain (for now).
        total_number_labels = size(Solver_setup.all_domain_labels,1);
        number_labels_per_domain = total_number_labels / Solver_setup.number_of_domains;
        Solver_setup.domain_labels = cell(Solver_setup.number_of_domains, number_labels_per_domain);

        % Keep track of the number of labels assigned per domain (should correspond with the 
        % above number)
        label_count_per_domain = zeros(Solver_setup.number_of_domains,1);

        % Loop now over all the domain labels and sort them correctly
        for label_index = 1:total_number_labels
            label = Solver_setup.all_domain_labels{label_index};
            % Extract the domain index portion, i.e. the @<ID> part, where <ID> is the domain id
            % to which this label belongs. We assume here we do not have @ in the name.
            domain_index = strsplit(label,{'@'});
            domain_index = str2num(domain_index{2});
            label_count_per_domain(domain_index) = label_count_per_domain(domain_index) + 1;
            if (label_count_per_domain(domain_index) > number_labels_per_domain)
                message_fc(Const,sprintf('Total number of labels per domain exceeded'));
                error(['Total number of labels per domain exceeded']);
            end%if
            Solver_setup.domain_labels(domain_index,label_count_per_domain(domain_index)) = cellstr(label);
        end%for
        
        % Allocate now space to store the basis functions list of each
        % domain (will be used when extracting e.g. matrix entries). We store two
        % types of vectors, one for keep track of the unknowns in the extended domain
        % (Solver_setup.rwg_basis_functions_domains) and the other for tracking only the 
        % internal RWGs (i.e. not included in the extended domain), 
        % Solver_setup.rwg_basis_functions_internal_domains
        Solver_setup.rwg_basis_functions_domains = cell(Solver_setup.number_of_domains,1);
        Solver_setup.rwg_basis_functions_internal_domains = cell(Solver_setup.number_of_domains,1);
        Solver_setup.rwg_basis_functions_unique_domains = cell(Solver_setup.number_of_domains,1);
        
        % Initialise now each of the domain basis function lists
        for domain_index = 1:Solver_setup.number_of_domains
            Solver_setup.rwg_basis_functions_domains{domain_index} = [];
            Solver_setup.rwg_basis_functions_unique_domains{domain_index} = [];
        end

        % Also keep track of all basis function allocation. Used below for unique allocation 
        % strategy. This is allocated below and initialsed to false.
        Solver_setup.rwg_basis_functions_allocated = zeros(Solver_setup.num_metallic_edges,1);
        
        % Let's assume we have disjoint domains (flag will be updated below
        % if we find this is not the case)
        Solver_setup.disconnected_domains = true;
    else
        % Only a single domain
        Solver_setup.number_of_domains = 1;
    end%if if (Const.domain_decomposition)
    
    % Initialise the return values:
    Solver_setup.rwg_basis_functions_shared_edge_nodes = zeros(Solver_setup.num_metallic_edges,2);    
    Solver_setup.rwg_basis_functions_trianglePlusFreeVertex = zeros(Solver_setup.num_metallic_edges,1);
    Solver_setup.rwg_basis_functions_triangleMinusFreeVertex = zeros(Solver_setup.num_metallic_edges,1);
    Solver_setup.rwg_basis_functions_shared_edge_centre = zeros(Solver_setup.num_metallic_edges,3);

    % Allocate an empty list to store the RWGs on the boundaries (just group all of them together)
    Solver_setup.rwg_basis_functions_on_interface = [];

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
        
            % Locate now the domain index for triangle plus and minus. Note, this is actually already 
            % given in the triangle label when we use the 'FA' card in FEKO to specify the domains @<ID>
            % see also comment above when we determine the unique set of labels that are associated with
            % each domain.
            % domain_index_triangle_plus  = find(strcmp(Solver_setup.domain_labels,triangle_plus_label));
            % domain_index_triangle_minus = find(strcmp(Solver_setup.domain_labels,triangle_minus_label));
            
            domain_index_triangle_plus = strsplit(triangle_plus_label,{'@'});
            domain_index_triangle_plus = str2num(domain_index_triangle_plus{2});
            domain_index_triangle_minus = strsplit(triangle_minus_label,{'@'});
            domain_index_triangle_minus = str2num(domain_index_triangle_minus{2});

            % If the basis function has not already been allocated to a domain, then do it now)
            if (~Solver_setup.rwg_basis_functions_allocated(edge_index))
                % Add it to the list of unique RWGs / domain
                Solver_setup.rwg_basis_functions_unique_domains{domain_index_triangle_plus} = [ ...
                    Solver_setup.rwg_basis_functions_unique_domains{domain_index_triangle_plus},...
                    edge_index];
                % Mark the RWG as allocated (i.e. assigned to a domain)
                Solver_setup.rwg_basis_functions_allocated(edge_index) = true;
            end%if

            % Add the basis function now to this domain's list
            if (domain_index_triangle_plus == domain_index_triangle_minus)
                % Basis function entirely in domain. Add it to the list.
                Solver_setup.rwg_basis_functions_domains{domain_index_triangle_plus} = [ ...
                    Solver_setup.rwg_basis_functions_domains{domain_index_triangle_plus},...
                    edge_index];                
            else
                % Basis function on the boundary - we have detected
                % interconnected domains!
                Solver_setup.disconnected_domains = false;
                
                % Add the basis function to both domains.
                Solver_setup.rwg_basis_functions_domains{domain_index_triangle_plus} = [ ...
                    Solver_setup.rwg_basis_functions_domains{domain_index_triangle_plus},...
                    edge_index];
                
                Solver_setup.rwg_basis_functions_domains{domain_index_triangle_minus} = [ ...
                    Solver_setup.rwg_basis_functions_domains{domain_index_triangle_minus},...
                    edge_index];

                Solver_setup.rwg_basis_functions_on_interface = [Solver_setup.rwg_basis_functions_on_interface edge_index];
            end%if

            if (LOCAL_DEBUG)
                fprintf(1,'Triangle %d (+) in domain %d\n',triangle_plus, domain_index_triangle_plus);
                fprintf(1,'Triangle %d (-) in domain %d\n',triangle_minus, domain_index_triangle_minus);
            end%if
            
        end%if (Const.domain_decomposition)
            
    end%for edge_index = 1:Solver_setup.num_metallic_edges

    % Loop over the domains and calculate the unknowns internal to each domain (i.e. excluding the RWG basis
    % functions on the interface). This is only necessary if we have interconnected domains:
    if (Const.domain_decomposition && ~Solver_setup.disconnected_domains)
        for el = 1:Solver_setup.num_finite_array_elements
            extended_domain_indices = Solver_setup.rwg_basis_functions_domains{el};    
            % Determine now this domain's internal unknowns, by just extracting 
            Solver_setup.rwg_basis_functions_internal_domains{el} = setdiff(extended_domain_indices, ...
                Solver_setup.rwg_basis_functions_on_interface);
        end%for
       
        % We need to extract now the domain interconnectivity here, so that
        % we can define the MBF primary generating sub-arrays later. Assume
        % for now a 1D case. TO-DO: Danie, we need to add an error check
        % here that makes sure we exit with an error should a 2D array be
        % detected. For 2D arrays as e.g. generated using the 'FA' card of FEKO, 
        % the domains are not always numbered 1,2,3,4,etc. for one row and then 
        % 5,6,7,8 for the next. We need to use e.g. the domain boundary information
        % to determine the interconnectivity.
        Solver_setup.domain_connectivity = 1:Solver_setup.num_finite_array_elements;

        % Set each of the domain types, i.e. from which generating sub-array to extract 
        % primary MBFs from (see also runMBFgenerator.m)
        %   Sub-array 1 : Left edge sub-array
        %   Sub-array 2 : Centre sub-array
        %   Sub-array 3 : Right edge sub-array
        Solver_setup.generating_subarray_parent = zeros(Solver_setup.num_finite_array_elements,1);
        Solver_setup.generating_subarray_parent(:)=2;   % Set all first to centre type
        Solver_setup.generating_subarray_parent(1) = 1; % Correct left edge
        Solver_setup.generating_subarray_parent(Solver_setup.num_finite_array_elements) = 3; % Correct right edge

        % Now that we have the interconnectivity, we define some information for the
        % generating sub-arrays that we will use to construct better primary MBFs
        if (Solver_setup.num_finite_array_elements > 3)

            % We have 3 generating sub-arrays for linear connected structures (which we consider at the moment)
            Solver_setup.generating_subarrays.number_of_domains = 3;            

            Solver_setup.generating_subarrays.domains = cell(3,1); % We will have 2 x edge domains and 1 centre domain
            Solver_setup.generating_subarrays.domains{1} = [Solver_setup.domain_connectivity(1), ...
                Solver_setup.domain_connectivity(2)];
            Solver_setup.generating_subarrays.domains{2} = [Solver_setup.domain_connectivity(2), ...
                Solver_setup.domain_connectivity(3), ...
                Solver_setup.domain_connectivity(4)];
            Solver_setup.generating_subarrays.domains{3} = [Solver_setup.domain_connectivity(Solver_setup.num_finite_array_elements-1), ...
                Solver_setup.domain_connectivity(Solver_setup.num_finite_array_elements)];
        else
            message_fc(Const,sprintf('Domain decomposition requires more than 3 array elements'));
            error(['Domain decomposition requires more than 3 array elements']);
        end

        % Extract now the edges associated with the sub-array domains:
        Solver_setup.generating_subarrays.rwg_basis_functions_domains = cell(3,1);
        for ii = 1:Solver_setup.generating_subarrays.number_of_domains
            subarray_domains = Solver_setup.generating_subarrays.domains{ii};            
            rwgs = [];
            % Loop over all the domains in this sub-array and extract the RWG elements.
            for jj = 1:length(subarray_domains)
                domain_index = subarray_domains(jj);
                rwgs = [rwgs Solver_setup.rwg_basis_functions_domains{domain_index}];
            end%for
            % Store the RWGs for this generating sub-array. Note: We
            % discard the RWGs on the interface between the domains, that
            % will be included in the Solver_setup.rwg_basis_functions_domains{jj} vector
            Solver_setup.generating_subarrays.rwg_basis_functions_domains{ii} = unique(rwgs);
        end%for

        % Also specify now the internal domain (in the domain_index list of the subarray) that will
        % the final (truncated) basis functions.
        Solver_setup.generating_subarrays.truncation_domain = cell(3,1); % We will have 2 x edge domains (left & right) and 1 centre domain
        
        domains = Solver_setup.generating_subarrays.domains{1};
        Solver_setup.generating_subarrays.truncation_domain{1} = domains(1); % Left edge sub-array (2 elements)
        
        domains = Solver_setup.generating_subarrays.domains{2};
        Solver_setup.generating_subarrays.truncation_domain{2} = domains(2); % Centre sub-array (3 elements)
        
        domains = Solver_setup.generating_subarrays.domains{3};
        Solver_setup.generating_subarrays.truncation_domain{3} = domains(2); % Right sub-array (2 elements)

    end %if
    