function [Z] = FillZMatrixByEdge(Const,Solver_setup)
    %function [Z] = FillZMatrixByEdge(omega,eps_0,mu_0,k,r_c,rho_c_pls,rho_c_mns,quad_pts,sing,dof2edge)
    %FillZMatrixByEdge
    %   Date: 2018.06.10
    %   Usage:
    %       [Z] = FillZMatrixByEdge(Const,Solver_setup)
    %
    %   Input Arguments:
    %       Const: 
    %           A global struct containing:
    %       Solver_setup
    %           A global struct containing solver setup details, e.g. frequency range,
    %           geometry data, basis function setup, etc.
    %   Output Arguments:
    %       Z
    %           The Z-matrices data calculated internally
    %
    %   Description:
    %       Fills the impedance matrix at each of the frequency samples.
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.253
    %   Stellenbosch University
    %   Email: dludick.sun.ac.za 
    %   Credit to: Prof. David B. Davidson for the FillZMatrixByEdge routine that is based on his
    %   MATLAB implementation as detailed in [1]
    %   
    %   References: 
    %   [1] David B. Davidson, Computational Electromagnetics for RF and Microwave Engineering, 
    %       Second Edition, (see Chapter 6)

    % Populate now the variables as used by [1] - esp. replace the use of
    % global variables.
    num_dofs = Solver_setup.num_metallic_edges;       % Replacing global NUM_DOFS
    elements = Solver_setup.triangle_vertices;        % Replacing global ELEMENTS
    node_coord = Solver_setup.nodes_xyz;              % Replacing global NODE_COORD
    ell = Solver_setup.rwg_basis_functions_length_m;  % Replacing global ELL
    
    quad_pts = Const.QUAD_PTS;
    sing     = Const.SING;
    eps_0    = Const.EPS_0;
    mu_0     = Const.MU_0;

    % Extract the triangle midpoints
    r_c = Solver_setup.triangle_centre_point;
    rho_c_pls = Solver_setup.rho_c_pls;
    rho_c_mns = Solver_setup.rho_c_mns;

    number_of_frequencies = Solver_setup.frequencies.freq_num; % Number of frequencies
    
    % Allocate some space for our impedance matrix
    Z = complex(zeros(num_dofs,num_dofs,number_of_frequencies)); % Generalized impedance matrix.

    % We will be calculating the Z matrix at each of the various frequencies:
    for freq_index = 1:number_of_frequencies

        % Extract the particular frequency value:
        freq = Solver_setup.frequencies.samples(freq_index);
        % Calculate some frequency dependent parameters required below
        omega = 2*pi*freq;  % Radial frequency
        lambda = Const.C0/freq;    % Wavelength in m
        k  = 2*pi/lambda;   % Wavenumber in rad/m

        message_fc(Const, sprintf('  Processing frequency %d of %d (%.2f Hz) ',freq_index,number_of_frequencies,freq))

        % Assemble by edges - not optimally fast computationally, but easy.
        % These are eqns. 32 and 33.
        for mm = 1:num_dofs

            for nn = 1:num_dofs
                % There are four terms to add here;
                % From integrals over source faces n+ & n- evaluated at centres of field faces m+ and m-.
                % The n integrals associate with q, the m, with p.
                % In datastructures EDGECONXELEMS & DOFLOCALNUM, second index 1
                % associates with +, 2 with -. 
                
                %pp_pls = EDGECONXELEMS(mm,1);            
                pp_pls = Solver_setup.rwg_basis_functions_trianglePlus(mm);
                %pp_mns = EDGECONXELEMS(mm,2);
                pp_mns = Solver_setup.rwg_basis_functions_triangleMinus(mm);

                %qq_pls = EDGECONXELEMS(nn,1);
                qq_pls = Solver_setup.rwg_basis_functions_trianglePlus(nn);
                %qq_mns = EDGECONXELEMS(nn,2);
                qq_mns = Solver_setup.rwg_basis_functions_triangleMinus(nn);
            
                % First, find contribution from n+ and n- faces evaluated at m+
                
                % --------------------------------------------------------------                
                % Look at Tm+
                % --------------------------------------------------------------
                % -- Contribution from Tn+
                triangle_tn_plus_free_vertex = Solver_setup.rwg_basis_functions_trianglePlusFreeVertex(qq_pls);
                [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell, pp_pls,qq_pls,nn,triangle_tn_plus_free_vertex,...
                    k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                Amn_pls_source_pls = MagVecPot;
                Phi_mn_pls_source_pls = -ScalPot;
                
                % -- Contribution from Tn-
                triangle_tn_minus_free_vertex = Solver_setup.rwg_basis_functions_triangleMinusFreeVertex(qq_mns);
                [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell, pp_pls,qq_mns,nn,triangle_tn_minus_free_vertex,...
                    k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                Amn_pls_source_mns = - MagVecPot;
                Phi_mn_pls_source_mns = +ScalPot;
                                
                Amn_pls = Amn_pls_source_pls + Amn_pls_source_mns;
                Phi_mn_pls = Phi_mn_pls_source_pls + Phi_mn_pls_source_mns;
                
                % --------------------------------------------------------------                
                % Look at Tm-
                % --------------------------------------------------------------
                % -- Contribution from Tn+
                [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell, pp_mns,qq_pls,nn,triangle_tn_plus_free_vertex,...
                    k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                Amn_mns_source_pls = MagVecPot;
                Phi_mn_mns_source_pls = -ScalPot;
                
                % -- Contribution from Tn-
                triangle_tn_minus_free_vertex = Solver_setup.rwg_basis_functions_triangleMinusFreeVertex(qq_mns);
                [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell, pp_mns,qq_mns,nn,triangle_tn_minus_free_vertex,...
                    k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                Amn_mns_source_mns = - MagVecPot;
                Phi_mn_mns_source_mns = +ScalPot;
                
                Amn_mns = Amn_mns_source_pls + Amn_mns_source_mns;
                Phi_mn_mns = Phi_mn_mns_source_pls + Phi_mn_mns_source_mns;
                       
                % Assemble with eq. 17 in [RWG82]
                Z(mm,nn) = 1i*omega*...
                    (dot(Amn_pls',rho_c_pls(mm,:))/2 + dot(Amn_mns',rho_c_mns(mm,:))/2) + Phi_mn_mns - Phi_mn_pls;
                
                %mm
                %nn
                Z(mm,nn) = ell(mm)* Z(mm,nn);
            end % for nn = 1:NUM_DOFS
        end %for mm = 1:NUM_DOFS
    end %for freq_index = 1:Solver_setup.frequencies.freq_num    
end %function FillZMatrixByEdge

% =================================================================================
% Local function for evaluating the potential integrals, as done in [1]
% =================================================================================
function [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell, field_pt,source_pt,source_edge, ii, k,r_c,quad_pts,sing,eps_0,mu_0,omega)
    % This subfunction computes the magnetic vector and scalar potentials for 
    % field point face field_pt and source point face source_pt. 

    % Code for debugging singularity scheme below:
    %[Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq_debug(field_pt,source_pt,r_c(field_pt,:),k,quad_pts,sing);
    [Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq(elements,node_coord, field_pt,source_pt,r_c(field_pt,:),k,quad_pts,sing);
    
    % Extract the nodes of the source triangle (q)
    qnodes = elements(source_pt,:);    
    r = zeros(3,3); % Store all three vertices of triangle q here in r()
    r(1,:) = node_coord(qnodes(1),:);
    r(2,:) = node_coord(qnodes(2),:);
    r(3,:) = node_coord(qnodes(3),:);
    
    % Extract the position of the ith free vertex (i.e. associated with the
    % ith RWG)
    %ii_nodes = elements(ii,:);
    rii = zeros(1,3);
    rii(1,1) = node_coord(ii,1);
    rii(1,2) = node_coord(ii,2);
    rii(1,3) = node_coord(ii,3);    
    
    %ii = DOFLOCALNUM(source_edge,source_tri); % This is the free vertex
    %associated with the source_edge - which is now passed here as an
    %argument by the calling routine.    
    MagVecPot = mu_0*ell(source_edge)/(4*pi)*...
        ( r(1,:)*Ipq_xi + r(2,:)*Ipq_eta + r(3,:)*Ipq_zeta - rii(1,:)*Ipq);
    ScalPot = ell(source_edge)/(1i*2*pi*omega*eps_0) * Ipq;
end % function Potentials
