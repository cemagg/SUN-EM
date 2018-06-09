function [dof_num,dof2edge] = renumber_RWG(dof_free_flag)
% RENUMBER_DOF renumbers the internal edges 
% (which are the degrees of freedom when using RWG elements).  
% The inverse mapping is also tracked - ie. which edge an RWG dof is 
% associated with.
% Author: D B Davidson, Dec 2009. Based on renumber_dof.
global NUM_EDGES NUM_DOFS
dof_num = zeros(1,NUM_EDGES); % 0 indicates prescribed dof.
counter = 0; 
for i_edge = 1:NUM_EDGES
   if (dof_free_flag(i_edge))
     counter = counter + 1;
     dof_num(i_edge) = counter;
     dof2edge(counter) = i_edge;
   end 
end
NUM_DOFS = counter;
