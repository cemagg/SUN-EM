function [x_nodes,y_nodes,z_nodes] = trimesh3D(a,b,Nx,Ny)
% TRIMESH3D Triangular mesher for plate, 3D coordinates. Produced a mesh with 2 Nx Ny
% elements. Assumes plate in z=0 plane.

global NODE_COORD ELEMENTS EDGES NUM_NODES NUM_ELEMENTS NUM_EDGES

NUM_ELEMENTS = 2 * Nx * Ny; 
% Following is the nodes the element consists of, numbered from lowest to
% highest nodes (used in subsequent EDGEMAKE function).
for jj = 1:Ny
    for ii = 1:Nx
        elem_offset = (jj-1)*2*Nx;
        node_offset = (jj-1)*(Nx+1);
        % L row.
        ELEMENTS(ii+elem_offset,:)    = [ii+node_offset    ii+1+node_offset     ii+Nx+1+node_offset];
        % "Inverted" L row.
        ELEMENTS(ii+elem_offset+Nx,:) = [ii+1+node_offset  ii+Nx+1+node_offset  ii+Nx+2+node_offset];
    end
end
   
NUM_NODES = (Nx+1)*(Ny+1);
% Following is nodal coordinates, (x y), per node.
del_x = a/Nx; 
del_y = b/Ny;
counter = 1;
for jj = 1:(Ny+1)
    for ii = 1:(Nx+1)
        %y_offset = (jj-1)*del_y
        NODE_COORD(counter,:) = [(ii-1)*del_x   (jj-1)*del_y   0];
        counter = counter+1;
    end
end

% for subsequent use by triplot.
x_nodes = [0 a/Nx a];
y_nodes = [0 b/Ny b];
z_nodes = zeros(1,NUM_NODES); % Redundant at present.
%keyboard;


