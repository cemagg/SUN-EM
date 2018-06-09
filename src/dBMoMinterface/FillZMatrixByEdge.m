function [Z] = FillZMatrixByEdge(omega,eps_0,mu_0,k,r_c,rho_c_pls,rho_c_mns,quad_pts,sing,dof2edge)
% FILLZMATRIXBYEDGE Fill the impedance matrix by edge pair.
% Code re-factored 24 Jan 2010, DBD. 

global ELEMENTS NODE_COORD NUM_DOFS EDGECONXELEMS DOFLOCALNUM ELL

Z = zeros(NUM_DOFS,NUM_DOFS); % Generalized impedance matrix.
% Assemble by edges - not optimally fast computationally, but easy.
% These are eqns. 32 and 33.
for mm = 1:NUM_DOFS
    for nn = 1:NUM_DOFS
        % There are four terms to add here;
        % From integrals over source faces n+ & n- evaluated at centres of field faces m+ and m-.
        % The n integrals associate with q, the m, with p.
        % In datastructures EDGECONXELEMS & DOFLOCALNUM, second index 1
        % associates with +, 2 with -. 
        
        pp_pls = EDGECONXELEMS(mm,1);
        pp_mns = EDGECONXELEMS(mm,2);
        qq_pls = EDGECONXELEMS(nn,1);
        qq_mns = EDGECONXELEMS(nn,2);
        
        % First, find contribution from n+ and n- faces evaluated at m+
        [MagVecPot,ScalPot] = Potentials(pp_pls,qq_pls,nn,1,k,r_c,dof2edge,quad_pts,sing,eps_0,mu_0,omega);
        Amn_pls_source_pls = MagVecPot;
        Phi_mn_pls_source_pls = -ScalPot;
        [MagVecPot,ScalPot] = Potentials(pp_pls,qq_mns,nn,2,k,r_c,dof2edge,quad_pts,sing,eps_0,mu_0,omega);
        Amn_pls_source_mns = - MagVecPot;
        Phi_mn_pls_source_mns = +ScalPot;
        Amn_pls = Amn_pls_source_pls + Amn_pls_source_mns;
        Phi_mn_pls = Phi_mn_pls_source_pls + Phi_mn_pls_source_mns;
        
        % Repeat now for m-
        [MagVecPot,ScalPot] = Potentials(pp_mns,qq_pls,nn,1,k,r_c,dof2edge,quad_pts,sing,eps_0,mu_0,omega);
        Amn_mns_source_pls = MagVecPot;
        Phi_mn_mns_source_pls = -ScalPot;
        [MagVecPot,ScalPot] = Potentials(pp_mns,qq_mns,nn,2,k,r_c,dof2edge,quad_pts,sing,eps_0,mu_0,omega);
        Amn_mns_source_mns = - MagVecPot;
        Phi_mn_mns_source_mns = +ScalPot;
        Amn_mns = Amn_mns_source_pls + Amn_mns_source_mns;
        Phi_mn_mns = Phi_mn_mns_source_pls + Phi_mn_mns_source_mns;
               
        % Assemble with eq. 17
        Z(mm,nn) = j*omega*...
            (dot(Amn_pls',rho_c_pls(mm,:))/2 + dot(Amn_mns',rho_c_mns(mm,:))/2) + Phi_mn_mns - Phi_mn_pls;
        
        %mm
        %nn
        Z(mm,nn) = ELL(dof2edge(mm))* Z(mm,nn);
    end
end

end

function [MagVecPot,ScalPot] = Potentials(field_pt,source_pt,source_edge,source_tri,k,r_c,dof2edge,quad_pts,sing,eps_0,mu_0,omega)
% This subfunction computes the magnetic vector and scalar potentials for 
% field point face field_pt and source point face source_pt. 
global ELEMENTS NODE_COORD DOFLOCALNUM ELL 

% Code for debugging singularity scheme below:
%[Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq_debug(field_pt,source_pt,r_c(field_pt,:),k,quad_pts,sing);
[Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq(field_pt,source_pt,r_c(field_pt,:),k,quad_pts,sing);
qnodes = ELEMENTS(source_pt,:);
r(1,:) = NODE_COORD(qnodes(1),:);
r(2,:) = NODE_COORD(qnodes(2),:);
r(3,:) = NODE_COORD(qnodes(3),:);
ii = DOFLOCALNUM(source_edge,source_tri);
MagVecPot = mu_0*ELL(dof2edge(source_edge))/(4*pi)*...
    ( r(1,:)*Ipq_xi + r(2,:)*Ipq_eta + r(3,:)*Ipq_zeta - r(ii,:)*Ipq);
ScalPot = ELL(dof2edge(source_edge))/(j*2*pi*omega*eps_0) * Ipq;
end
