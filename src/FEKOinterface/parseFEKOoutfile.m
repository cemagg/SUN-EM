function [Const, FEKO_data] = parseFEKOoutfile(Const)
    %parseFEKOoutfile
    %   Date: 2017.07.01
    %   Usage:
    %       [Const, FEKO_data] = parseFEKOoutfile(Const)
    %
    %   Input Arguments:
    %       Const: A global struct containing:
    %       FEKOoutfilename
    %           FEKO *.out filename (e.g. 'yagi.out')
    %
    %   Output Arguments:
    %       FEKO_data:
    %           The struct containing the frequency list data (e.g. number
    %           of points, frequency sample values, etc.), triangle information and
    %           basis function setup information.
    %
    %   Description:
    %       Reads in a FEKO *.out file and extracts the FEKO data (frequencies, triangle info, 
    %       basis function setup).
    %
    %   =======================
    %   Written by Danie Ludick on June 13, 2013
    %   Stellenbosch University
    %   Email: dludick.sun.ac.za

    %   Please note that additional information on reading the *.out files using FEKO can be found 
    %   on the FEKO website:
    %   http://www.feko.info/
    
    narginchk(1,1);

    fid = fopen(Const.FEKOoutfilename,'r');

    if fid == -1
        message_fc(Const,sprintf('Error reading FEKO *.out file: %s',Const.FEKOoutfilename));
        error(['Error reading FEKO *.out file: %s' Const.FEKOoutfilename]);
    end

    message_fc(Const,' ');
    message_fc(Const,...
        '------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Reading additional FEKO data from *.out file'));
    message_fc(Const,sprintf('  *.out file: %s',Const.FEKOoutfilename) );

    % ================================
    % Initialisations
    % ================================
    FEKO_data = [];
    FEKO_data.frequencies = [];    
    %FEKO_data.rwg_basis_functions = [];
    
    end_flag = 0;
    FEKO_data.frequencies.samples = [];
    freq_num = 0;
    FEKO_data.num_finite_array_elements = -1;

    % A local debug flag (e.g. to plot geometry)
    LOCAL_DEBUG = false;

    while end_flag == 0
        line=fgetl(fid);
        
        % Check for end of file:
        if strcmp(line,'                    SUMMARY OF REQUIRED TIMES IN SECONDS')
            end_flag = 1;
        end%if
        
        % -------------------------------------------------
        % -- Extract the frequency sample data
        % Note: for multiple solutions, this information will be written to the *.out file
        % for each. We therefore have to make sure that only the unique freq. increments are
        % retained - this is done after the file has been parsed.
        g = strfind(line,'Frequency in Hz');        
        if (g > 0)
            freq_num = freq_num + 1;
            %line
            freq_inc_info = strsplit(line);

            % Loop now over each of the elements in this frequency increment cell
            % array and extract the numeric value for the frequency data
            for cell_data = 1:size(freq_inc_info,2)
                str = freq_inc_info{cell_data};
                if (~isempty(str2num(str)))
                    FEKO_data.frequencies.samples(freq_num) = str2num(str);
                end%if

            end%for
        end%if (frequency increment)

        % -------------------------------------------------
        % -- Parse the FR card
        g = strfind(line,'FR ');
        if (g > 0)
            freq_card_info = strsplit(line); % 'CARD NAME',i1,i2,i3,i4,i5,f1,f2,f3,f4,f5
            num_freq_samples = str2num(freq_card_info{2});
            freq_start = str2num(freq_card_info{7});
            freq_end = str2num(freq_card_info{9});
        end% if 'FR' card

        % -------------------------------------------------
        % -- Parse the number of metallic triangles
        g = strfind(line,'Number of metallic triangles:'); 
        if (g > 0)
            triangle_line_info = strsplit(line);
            FEKO_data.num_metallic_triangles = str2num(triangle_line_info{6});
        end%if

        % -------------------------------------------------
        % -- Parse the number of metallic edges
        g = strfind(line,'Number of metallic edges (MoM):'); 
        if (g > 0)
            edge_line_info = strsplit(line);
            FEKO_data.num_metallic_edges = str2num(edge_line_info{7});
        end%if

        % -------------------------------------------------
        % -- Parse the number of array elements 
        % Note: Hopefully this is the only reference of this sort
        g = strfind(line,'Number of elements'); 
        if (g > 0)
            array_line_info = strsplit(line);
            FEKO_data.num_finite_array_elements = str2num(array_line_info{6});
        end%if

    end%while end_flag == 0
    fclose('all');

    % ==============================================
    % Some additional processing for frequency data
    % ==============================================
    % Correct now the frequency samples to take into account that multiple solutions / frequency
    % might have been specified (see above comment also)
    FEKO_data.frequencies.samples = unique(FEKO_data.frequencies.samples);
    FEKO_data.frequencies.freq_num = size(FEKO_data.frequencies.samples,2);

    % The following variables are used in postCMA
    Const.freqStart = freq_start;
    Const.freqEnd   = freq_end;
    Const.freqData  = FEKO_data.frequencies.samples;
    Const.numFreq   = length(Const.freqData);

    % Check here that the number of frequency samples read, corresponds to the number specified
    % in the 'FR' card
    if (FEKO_data.frequencies.freq_num ~= num_freq_samples)
        message_fc(Const,sprintf('[parseFEKOoutfile]: Number of frequencies read incorrect'));
        error(['[parseFEKOoutfile]: Number of frequencies read incorrect']);
    end%if
    
    % ==============================================
    % Some additional processing for triangle data
    % ==============================================
    % Now that we know how many metallic triangles we have, we can loop over this block in the 
    % *.out file and extract that information

    % Open the *,out file again
    fid = fopen(Const.FEKOoutfilename,'r');

    if fid == -1
        message_fc(Const,sprintf('Error reading FEKO *.out file: %s',Const.FEKOoutfilename));
        error(['Error reading FEKO *.out file: %s' Const.FEKOoutfilename]);
    end

    end_flag = 0;

    % Create an empty cell array to store the (variable) string labels of each triangle
    FEKO_data.metallic_triangles_labels = cell(FEKO_data.num_metallic_triangles,1);

    % Create other datastructures to store the triangle information. Note we initialise
    % this as NaN, otherwise we might run into a problem if the first node is also
    % at position (0,0,0). Note, the node_xyz matrix will also contain node indices
    % of other structures later, e.g. that of segments, or triangle + segment interfaces.
    FEKO_data.nodes_xyz = nan(FEKO_data.num_metallic_triangles,3);

    % -- Store only the nodal indices here below. Actual cartesian co-ordinates of the nodes
    % are stored in the nodes_xyz structure.
    FEKO_data.triangle_vertices = zeros(FEKO_data.num_metallic_triangles,3);

    FEKO_data.triangle_area_m2 = zeros(FEKO_data.num_metallic_triangles,1);

    % -- Initialisations for the RWG basis functions (Type 1 elements as per FEKO convention)
    FEKO_data.rwg_basis_functions_length_m = zeros(FEKO_data.num_metallic_edges,1);
    FEKO_data.rwg_basis_functions_triangleMinus = zeros(FEKO_data.num_metallic_edges,1);
    FEKO_data.rwg_basis_functions_trianglePlus = zeros(FEKO_data.num_metallic_edges,1);
    
    while end_flag == 0
        line=fgetl(fid);
        
        % Check for end of file:
        if strcmp(line,'                    SUMMARY OF REQUIRED TIMES IN SECONDS')
            end_flag = 1;
        end%if
        
        % -------------------------------------------------
        % -- Extract the triangle data
        g = strfind(line,'DATA OF THE METALLIC TRIANGLE');
        if (g > 0)
            % 6 lines before start of triangle data
            for ii=1:6
                line=fgetl(fid);    
            end%for
            
            max_node_index = 0;
            % Read now the data for all the triangles
            for ii=1:FEKO_data.num_metallic_triangles
                %fprintf('    ** processing triangle %d\n', ii);
                %line
                triangle_line_data = strsplit(line);
                triangle_id = str2num(triangle_line_data{2});
                FEKO_data.metallic_triangles_labels{triangle_id} = triangle_line_data{3};
                
                % 2018.05.15: Read now the other triangle data, e.g. nodal co-ordinates
                % area, shared edges info.

                % Read the X, Y and Z nodal co-ordinates of the triangle
                % Note: might still contain duplicate entries at this point.

                % Read 3 nodes:
                for jj=1:3                   
                    node_xyz = [str2num(triangle_line_data{4}), ...
                    str2num(triangle_line_data{5}), ...
                    str2num(triangle_line_data{6})];

                    % Take care before inserting this node address into the nodal list.
                    % We only store unique node identifiers, not duplicate nodes. We 
                    % check therefore first whether the node is in the list. If it is, 
                    % then the routine below returns that index.

                    node_index = is_node_in_list(Const, node_xyz, FEKO_data.nodes_xyz);

                    % We have detected a new node - add this to the back of the list.
                    if ( node_index == -1)
                        % New node - add to back of list
                        max_node_index = max_node_index + 1;
                        node_index = max_node_index;
                    end%if

                    % Update the node and triangle data structures.
                    FEKO_data.nodes_xyz(node_index,:) = node_xyz;
                    FEKO_data.triangle_vertices(triangle_id,jj) = node_index;

                    % Read the next node
                    line=fgetl(fid);
                    triangle_line_data = strsplit(line);

                    % TO-DO: Later, we can also extract the medium on both sides of the
                    % triangle. This is associated with the second and third nodal co-ordinate entry.

                end%for
                
                % Extract the normal co-ordinate of the triangle
                FEKO_data.triangle_normal_vector(triangle_id,:) = [str2double(triangle_line_data{2}), ...
                    str2double(triangle_line_data{3}), ...
                    str2double(triangle_line_data{4})];

                % Extract the triangle area (in m^2)
                FEKO_data.triangle_area_m2(triangle_id,1) = str2double(triangle_line_data{5});

                % Extract the triangle centre-point [X,Y,Z] co-ordinate by first extracting the 
                % triangle vertices and then calculating from that the midpoint.
                v1 = FEKO_data.triangle_vertices(triangle_id,1);
                v2 = FEKO_data.triangle_vertices(triangle_id,2);
                v3 = FEKO_data.triangle_vertices(triangle_id,3);
    
                v1XYZ = [FEKO_data.nodes_xyz(v1,1), FEKO_data.nodes_xyz(v1,2), FEKO_data.nodes_xyz(v1,3)];
                v2XYZ = [FEKO_data.nodes_xyz(v2,1), FEKO_data.nodes_xyz(v2,2), FEKO_data.nodes_xyz(v2,3)];
                v3XYZ = [FEKO_data.nodes_xyz(v3,1), FEKO_data.nodes_xyz(v3,2), FEKO_data.nodes_xyz(v3,3)];

                FEKO_data.triangle_centre_point(triangle_id,:) = (1/3) .* (v1XYZ + v2XYZ + v3XYZ);

                % Read the next triangle entry
                line=fgetl(fid);
            end

        end% if g>0

        % -------------------------------------------------
        % -- Extract the triangle edge data (RWG basis functions)
        % -------------------------------------------------        
        g = strfind(line,'DATA OF THE METALLIC EDGES');
        if (g > 0)
            % 5 lines before start of triangle data
            for ii=1:4
                line=fgetl(fid);
            end%for

            % Read now the data for all the metallic edges
            for ii=1:FEKO_data.num_metallic_edges
                %fprintf('    ** processing RWG basis function %d\n', ii);
                %line                

                rwg_bf_data = strsplit(line);         

                rwg_bf_id   = str2num(rwg_bf_data{2});
                % Check that we do not exceed the max. number of BFs.
                if (rwg_bf_id > FEKO_data.num_metallic_edges)
                    message_fc(Const,sprintf('Number of RWG basis functions exceeded'));
                    error(['Number of RWG basis functions exceeded']);
                end%if

                rwg_bf_type = str2num(rwg_bf_data{3});
                % At the moment, we only allow Type 1 BFs, i.e. RWG elements spanning
                % the shared edge between two triangles
                if (rwg_bf_type > 1)
                    message_fc(Const,sprintf('Only RWG basis functions supported'));
                    error(['Only RWG basis functions supported']);
                end%if

                % Extract the edge length [m]
                FEKO_data.rwg_basis_functions_length_m(rwg_bf_id) = str2double(rwg_bf_data{4});

                % Extract the Tm+ triangle (KORP in FEKO) [m]
                FEKO_data.rwg_basis_functions_trianglePlus(rwg_bf_id) = str2num(rwg_bf_data{8});

                % Extract the Tm- triangle (KORM in FEKO) [m]
                FEKO_data.rwg_basis_functions_triangleMinus(rwg_bf_id) = str2num(rwg_bf_data{9});

                % Read the next basis function entry
                line=fgetl(fid);
            end% for ii=1:FEKO_data.num_metallic_edges
        end% if g>0

    end%while end_flag == 0
    fclose('all');

    % ==============================================
    % Additional processing for the basis functions
    % ==============================================

    % Extract the free vertex of Tm+ and Tm- for the RWGs
    % and the nodal indices of the shared edge of the RWG
    FEKO_data = extract_shared_edge_triangle_details(Const, FEKO_data);

    % Set the total number of MoM basis functions
    % -- For now, we have RWG elements. Add as more types are included.
    FEKO_data.num_mom_basis_functions = FEKO_data.num_metallic_edges;

    % ==============================================
    % Additional processing for the finite antenna array (if present)
    % ==============================================
    if (FEKO_data.num_finite_array_elements ~= -1)
        % TO-DO: Note, at this stage we assume that the entire MoM subset is
        % allocated to the finite array.
        FEKO_data.mom_basis_functions_per_array_element = ...
            FEKO_data.num_mom_basis_functions/FEKO_data.num_finite_array_elements;
    else
        FEKO_data.mom_basis_functions_per_array_element = -1;
    end%if

    % ==============================================
    % Some local debugging - display the mesh
    % ==============================================
    if (LOCAL_DEBUG)
        displayTriangleMesh(Const, FEKO_data.triangle_vertices, FEKO_data.triangle_centre_point, FEKO_data.nodes_xyz);
    end%if

    % ==========================================================================================
    % DONE : Print a summary
    message_fc(Const,sprintf('\n  Number of frequencies: %d',FEKO_data.frequencies.freq_num));    
    message_fc(Const,sprintf('  Number of metallic triangles: %d',FEKO_data.num_metallic_triangles));
    message_fc(Const,sprintf('  Number of metallic edges (RWG basis functions): %d',FEKO_data.num_metallic_edges));
    message_fc(Const,sprintf('  Number of MoM basis functions: %d',FEKO_data.num_mom_basis_functions));
    message_fc(Const,sprintf('  Number of finite array elements: %d',FEKO_data.num_finite_array_elements));
    message_fc(Const,sprintf('  Number of MoM basis functions per array element : %d', ...
        FEKO_data.mom_basis_functions_per_array_element));
    message_fc(Const,sprintf('  Finished processing the *.out file: %s',Const.FEKOoutfilename));