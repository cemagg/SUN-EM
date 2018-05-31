function node_index = is_node_in_list(Const, node_xyz, node_list_xyz)
    %lineCentre
    %   Usage:
    %       node_index = is_node_in_list(node_xyz, node_list)
    %
    %   Input Arguments:
    %       Const
    %           Constant containing some parameters, e.g. eps
    %       node_xyz
    %           The cartesian co-ordinates of a node
    %       node_list
    %           The cartesian co-ordinates of a list of nodes
    %   Output Arguments:
    %       node_index
    %           Index of the node_xyz located in node_list_xyz, 
    %           otherwise -1 if it cannot be found.
    %   Description:
    %       Checks whether a particular nodal co-ordinate (XYZ) 
    %       is located in a list of nodes. If so, then the index
    %       in the list is returned. If not, then -1 is returned.
    %   =======================
    %   Written by Danie Ludick on 2018.05.15
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    error(nargchk(3,3,nargin));

    % Initialisation
    node_index = -1;

    % loop over the node_list_xyz and check each X, Y and Z co-ordinate
    for ii=1:length(node_list_xyz)
        if ( (abs(node_xyz(1) - node_list_xyz(ii,1)) <= Const.EPS) && ...
             (abs(node_xyz(2) - node_list_xyz(ii,2)) <= Const.EPS) && ...
             (abs(node_xyz(3) - node_list_xyz(ii,3)) <= Const.EPS) )
            % -- The node has been found
            node_index = ii;
            return;
        end%if
    end%for