function edgemake_MoM
% EDGEMAKE_MoM Make edges on triangular mesh.
% Based on the 2D FEM code.
% DB Davidson, Dec 2009.


global ELEMENTS ELEMENT_EDGES NUM_ELEMENTS EDGES NUM_EDGES LOCALEDGENODES ELL ...
    NODE_COORD 
% Element 1 first:

EDGES(1,:) = [ELEMENTS(1,LOCALEDGENODES(1,1)),ELEMENTS(1,LOCALEDGENODES(1,2))];
EDGES(2,:) = [ELEMENTS(1,LOCALEDGENODES(2,1)),ELEMENTS(1,LOCALEDGENODES(2,2))];
EDGES(3,:) = [ELEMENTS(1,LOCALEDGENODES(3,1)),ELEMENTS(1,LOCALEDGENODES(3,2))];
ELEMENT_EDGES(1,1) = 1;
local_edges(1,1)  = 1;
ELEMENT_EDGES(1,2) = 2;
local_edges(2,1)  = 2;
ELEMENT_EDGES(1,3) = 3;
local_edges(3,1)  = 3;
% Now other elements
edge_counter = 3;
for ielem = 2:NUM_ELEMENTS
   for jedge = 1:3
     TEMPEDGES(1,:) = [ELEMENTS(ielem,LOCALEDGENODES(jedge,1)),ELEMENTS(ielem,LOCALEDGENODES(jedge,2))];
     new_edge = 1; % Default: true. Re-set for each edge.
     % Now test if this edge has already been assigned
     for kedge = 1:edge_counter
       if (TEMPEDGES(1,:) == EDGES(kedge,:)) 
         new_edge = 0;
         ELEMENT_EDGES(ielem,jedge) = kedge;
         break;
       end
     end
     if new_edge 
       edge_counter = edge_counter+1;
       EDGES(edge_counter,:) = TEMPEDGES(1,:);
       ELEMENT_EDGES(ielem,jedge) = edge_counter;
     end
   end
end
NUM_EDGES = edge_counter;
% Find edge lengths
for iedge = 1:NUM_EDGES
    ELL(iedge) = norm(NODE_COORD(EDGES(iedge,2),:)-NODE_COORD(EDGES(iedge,1),:));
end
