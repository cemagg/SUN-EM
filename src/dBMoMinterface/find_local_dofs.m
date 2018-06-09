function find_local_dofs(dof_RWG)
%FIND_LOCAL_DOFS resolves the local edge numbers for the two faces shared by each edge-based dof and
%saves them in datastructure DOFLOCALNUM
% DB Davidson, Dec 09

global  ELEMENT_EDGES NUM_ELEMENTS NUM_EDGES DOFLOCALNUM NUM_DOFS

DOFLOCALNUM  = zeros(NUM_DOFS,2);

for iedge = 1:NUM_EDGES
    edge_found = 0; % flag
    terminate = 0;
    for jelem = 1:NUM_ELEMENTS
        for kedge = 1:3
            if ELEMENT_EDGES(jelem,kedge) == iedge
                if ~edge_found
                    if dof_RWG(iedge) % i.e. this edge has a dof associated with it
                        DOFLOCALNUM(dof_RWG(iedge),1) = kedge;
                    end
                    edge_found = 1;
                else
                    if dof_RWG(iedge)
                        DOFLOCALNUM(dof_RWG(iedge),2) = kedge;
                    end
                    terminate = 1; % Terminate element loop; max of two faces per edge
                    break
                end
            end
        end
        if terminate
            break % exit element loop and return to outermost edge loop.
        end
    end
end

end

