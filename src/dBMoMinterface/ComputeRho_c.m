function [rho_c_pls,rho_c_mns,x_grid_pls,y_grid_pls,x_grid_mns,y_grid_mns] = ComputeRho_c(r_c)
% COMPUTERHO_C computes the vectors \vec\rho_n^c+ and - in [Fig2, RWG82].

global ELEMENTS NODE_COORD NUM_DOFS EDGECONXELEMS DOFLOCALNUM LOCALVERTEX

rho_c_pls = zeros(NUM_DOFS,3); 
rho_c_mns = zeros(NUM_DOFS,3); 

x_grid_pls = zeros(NUM_DOFS,1);
y_grid_pls = zeros(NUM_DOFS,1);
for mm=1:NUM_DOFS % General code
    pp_pls = EDGECONXELEMS(mm,1);
    loc_node_pls = LOCALVERTEX(DOFLOCALNUM(mm,1));
    vertx = NODE_COORD(ELEMENTS(pp_pls,loc_node_pls),:);
    %r_c(pp_pls,:)
    rho_c_pls(mm,:) = r_c(pp_pls,:) - vertx; % Directed from vertex
    x_grid_pls(mm) = r_c(pp_pls,1);
    y_grid_pls(mm) = r_c(pp_pls,2);
    pp_mns = EDGECONXELEMS(mm,2);
    loc_node_mns = LOCALVERTEX(DOFLOCALNUM(mm,2));
    vertx = NODE_COORD(ELEMENTS(pp_mns,loc_node_mns),:);
    rho_c_mns(mm,:) = - (r_c(pp_mns,:) - vertx); % Directed to vertex
    x_grid_mns(mm) = r_c(pp_mns,1);
    y_grid_mns(mm) = r_c(pp_mns,2);
    
    % Note that on the triangular from rectangular mesh, rho_c_pls and 
    % rho_c_mns are co-linear, but this is not generally true.
end


end

