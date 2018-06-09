function edge_conx_elem(dof_RWG)
% EDGE_CONX_elem sets up a connection matrix, which lists the two elements
% each interior edge (the degree of freedom) is connected to.
% EDGECONXELEMS(i,j) indicates that the dof associated with edge i is
% connected to element j.
% DB Davidson, Dec 09.

global EDGECONXELEMS ELEMENT_PLS_MNS NUM_EDGES NUM_ELEMENTS ELEMENT_EDGES NUM_DOFS

% An interior edge must be shared by two triangles
EDGECONXELEMS = zeros(NUM_DOFS,2);
for iedge = 1:NUM_EDGES
    if dof_RWG(iedge) % this is an interior edge
        counter = 1;
        for jelem = 1:NUM_ELEMENTS
            for kedge = 1:3
                if ELEMENT_EDGES(jelem,kedge) == iedge
                    EDGECONXELEMS(dof_RWG(iedge),counter) = jelem;
                    if counter == 1
                        ELEMENT_PLS_MNS(jelem,kedge) = +1;
                    else
                        ELEMENT_PLS_MNS(jelem,kedge) = -1;
                    end
                    counter = counter+1;
                end
            end
        end
    end
end
