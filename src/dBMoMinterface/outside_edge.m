function dof_flag = outside_edge(a,b)
% OUTSIDE_EDGE flags the edges as internal (1) or external (0). 
% Author: D B Davidson, December 2009.
% Based on FEM code dof_free.
% Note that this code assumes the scatterer is a rectangular plate, lying a
% plane of constant z, with one corner on the origin and the other at
% (x=a,y=b).

global NUM_EDGES EDGES NODE_COORD
dof_flag(1:NUM_EDGES) = 1; % Flag all as internal to start.

for i_edge = 1:NUM_EDGES
   node1 = EDGES(i_edge,1);
   node2 = EDGES(i_edge,2);
   if abs(NODE_COORD(node1,2))<eps  && abs(NODE_COORD(node2,2)) <eps % ie y=0
       dof_flag(i_edge) = 0;
   end
   if abs(NODE_COORD(node1,2)-b)<eps  && abs(NODE_COORD(node2,2)-b) <eps % ie y=b
       dof_flag(i_edge) = 0;
   end
   if abs(NODE_COORD(node1,1))<eps  && abs(NODE_COORD(node2,1)) <eps % ie x=0
       dof_flag(i_edge) = 0;
   end
   if abs(NODE_COORD(node1,1)-a)<eps  && abs(NODE_COORD(node2,1)-a) <eps % ie x=a
       dof_flag(i_edge) = 0;
   end
end