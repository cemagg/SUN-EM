function [Z] = FillZMatrixByFace(omega,eps_0,mu_0,k,r_c,rho_c_pls,rho_c_mns,quad_pts,sing,dof2edge,dof_RWG)
% FILLZMATRIXBYFACE Fill the impedance matrix by face pair.
%  Code computes the integrals for one source triangle-field triangle
%  interaction and then assembles by element, similar to a FEM code.
%  Filling by face is somewhat faster than by edge, but not the factor of
%  nine suggested in [RWG82]. The reason is probably that the code is more
%  complex. Gains are also more noticeable with more quadrature points,
%  where the overheads of this scheme impact less.


global ELEMENTS NODE_COORD NUM_DOFS EDGECONXELEMS DOFLOCALNUM ELL NUM_ELEMENTS ...
    ELEMENT_EDGES LOCALVERTEX ELEMENT_PLS_MNS

Z = zeros(NUM_DOFS,NUM_DOFS); % Generalized impedance matrix.

% There are four terms to add here;
% From integrals over source faces n+ & n- evaluated at centres of field faces m+ and m-.
% The n integrals associate with q, the m, with p.
% In datastructures EDGECONXELEMS & DOFLOCALNUM, second index 1
% associates with +, 2 with -.
% This must be done for each of the three interior edges per element.

% These are eqns. 32 and 33.
for melem = 1:NUM_ELEMENTS
    for nelem = 1:NUM_ELEMENTS
        % Compute the potentials (for the magnetic vector potential, for all three nodes)
        % - only computed once per combination of source and field faces.
        [MagVecPot,ScalPot] = Potentials(melem,nelem,k,r_c,quad_pts,sing,eps_0,mu_0,omega);
        % Scale by source edge
        for nedge = 1:3
            nn = dof_RWG(ELEMENT_EDGES(nelem,nedge)); % Find the corresponding dof for this source edge
            if nn % Only compute if source point is an interior edge
                nnode = LOCALVERTEX(nedge);
                % Note that distinguishing between node and edge is not
                % necessary with current local vertex numbering scheme
                % (which is an identity mapping), but is included for
                % completeness.
                EllMagVecPot(:,nnode) = ELL(dof2edge(nn))*MagVecPot(:,nnode);
                EllScalPot   = ELL(dof2edge(nn))*ScalPot;
                % Now resolve whether this is a + or - source triangle for this edge dof
                % (note that a triangle can be + for one edge, but - for another).
                Amn = ELEMENT_PLS_MNS(nelem,nedge)*EllMagVecPot(:,nnode);
                Phi_mn = -ELEMENT_PLS_MNS(nelem,nedge)*EllScalPot;
                % Assemble into eq(17),[RGW82]
                for medge = 1:3
                    mm = dof_RWG(ELEMENT_EDGES(melem,medge)); % Find the corresponding dof for this field edge
                    if mm % only compute if field point is an interior edge
                        if ELEMENT_PLS_MNS(melem,medge) == +1
                            rho_c = rho_c_pls(mm,:);
                        else
                            rho_c = rho_c_mns(mm,:);
                        end
                        Z(mm,nn) = Z(mm,nn)+ ELL(dof2edge(mm))*(j*omega*...
                            dot(Amn',rho_c)/2 - Phi_mn*ELEMENT_PLS_MNS(melem,medge));
                    end
                end
            end
        end
    end
    
end

end
function [MagVecPot,ScalPot] = Potentials(field_pt,source_pt,k,r_c,quad_pts,sing,eps_0,mu_0,omega)
% This subfunction computes the magnetic vector and scalar potentials for
% field point face field_pt and source point face source_pt. It returns the
% former for all three possible integrals in [eq.32, RWG82]
global ELEMENTS NODE_COORD DOFLOCALNUM ELL

[Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq(field_pt,source_pt,r_c(field_pt,:),k,quad_pts,sing);
qnodes = ELEMENTS(source_pt,:);
r(1,:) = NODE_COORD(qnodes(1),:);
r(2,:) = NODE_COORD(qnodes(2),:);
r(3,:) = NODE_COORD(qnodes(3),:);
for node = 1:3
    MagVecPot(:,node) = mu_0/(4*pi)*...
        ( r(1,:)*Ipq_xi + r(2,:)*Ipq_eta + r(3,:)*Ipq_zeta - r(node,:)*Ipq);
end
ScalPot = 1/(j*2*pi*omega*eps_0) * Ipq;
end
