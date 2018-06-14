function [Z] = FillZMatrixByEdge(Const,Solver_setup)    
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
    %       Second Edition, (see Chapter 6) - [DBD2011]
    %   [2] Xinlei Chen, Changqing Gu, Zhenyi Niu, and Zhuo L, "Fast Dipole Method for Electromagnetic Scattering From Perfect 
    %       Electric Conducting Targets", IEEE TRANSACTIONS ON ANTENNAS AND PROPAGATION, VOL. 60, NO. 2, FEBRUARY 2012

    message_fc(Const, sprintf('  Calculating Z-matrix '));
    
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

    % Set some general parameters
    number_of_frequencies = Solver_setup.frequencies.freq_num; % Number of frequencies
    Z.numFreq = number_of_frequencies;
    Z.mBasis  = num_dofs; % number of rows (testing/field functions)
    Z.nBasis  = num_dofs; % number of cols (basis/source  functions)
    
    % Allocate some space for our impedance matrix
    Z.values = complex(zeros(num_dofs,num_dofs,number_of_frequencies)); % Generalized impedance matrix.

    % We will be calculating the Z matrix at each of the various frequencies:
    for freq_index = 1:number_of_frequencies

        % Extract the particular frequency value:
        freq = Solver_setup.frequencies.samples(freq_index);
        % Calculate some frequency dependent parameters required below
        omega = 2*pi*freq;       % Radial frequency
        lambda = Const.C0/freq;  % Wavelength in m
        k  = 2*pi/lambda;        % Wavenumber in rad/m
        jk = 1i*k;

        if (Const.useEDM)
            % Precompute a constant for the Equivalent Dipole Method (EDM) (see [2])
            C1 = Const.ETA_0/(4*pi);
            threshold_EDM = 1.5*lambda;  % Set the EDM threshold. In e.g. [2], this is selected as 0.15*lambda
        end%if (Const.useEDM)

        message_fc(Const, sprintf('    Processing frequency %d of %d (%.2f Hz) ',freq_index,number_of_frequencies,freq))

        % Assemble by edges - not optimally fast computationally, but easy.
        % These are eqns. 32 and 33.
        num_edm_approximations = 0;
        num_std_zmn_entries_calculated = 0;
        for mm = 1:num_dofs

            %pp_pls = EDGECONXELEMS(mm,1);            
            pp_pls = Solver_setup.rwg_basis_functions_trianglePlus(mm);
            %pp_mns = EDGECONXELEMS(mm,2);
            pp_mns = Solver_setup.rwg_basis_functions_triangleMinus(mm);

            if (Const.useEDM)
                % Extract the precalculated dipole moment for edge mm
                ed_m_moment = Solver_setup.rwg_basis_functions_equivalent_dipole_moment(mm);
                % Centre of the above equivalent dipole
                r_ed_m_centre = Solver_setup.rwg_basis_functions_equivalent_dipole_centre(mm);
            end
           
            for nn = 1:num_dofs

                % 2018.06.13: We can use the Equivalent Dipole Method (EDM), if active and only if the distance 
                % between the two RWG dipole moments is sufficent (otherwise we still use normal matrix fill below)
                Zmn_calculated = false;
                if (Const.useEDM)

                    % Extract the precalculated dipole moment for edge mm
                    ed_n_moment = Solver_setup.rwg_basis_functions_equivalent_dipole_moment(nn);
                    % Centre of the above equivalent dipole
                    r_ed_n_centre = Solver_setup.rwg_basis_functions_equivalent_dipole_centre(nn);

                    % Check first whether the distance between the dipole moments are sufficiently large (i.e. larger
                    % than the preset threshold above). If it is, then we can use the EDM - if not, stil perform standard
                    % matrix calculation below.
                    R = r_ed_m_centre - r_ed_n_centre;
                    R_mag = abs(R);
                    one_over_R_mag = 1/R_mag;
                    R_hat = R*one_over_R_mag; % normalised R vector

                    if ((R_mag > threshold_EDM))                        
                        num_edm_approximations = num_edm_approximations + 1;

                        % Calculate C - see [2, Eq. 7]. Trying to do it a bit more efficiently by using precalculated terms
                        % such as one_over_R_mag. TO-DO: Test the impact of this ... vs. readability of code.
                        C = (one_over_R_mag)*(one_over_R_mag)*(1 + (1/jk)*one_over_R_mag);

                        % Calculate Zmn as done in [2, Eq. 6]
                        Z.values(mm,nn) = C1*exp(-jk*R_mag) * ( dot(ed_m_moment,ed_n_moment)*(jk*one_over_R_mag + C) - ...
                             dot(ed_m_moment,R_hat)*dot(R_hat,ed_n_moment)*(jk*one_over_R_mag + 3*C) );

                        % Set a flag that we are now done with this matrix entry
                        Zmn_calculated = true; % using EDM
                    end%if

                end % if (Const.useEDM)
                    
                % Perform the standard Matrix calculation (using Gausian quadrature integration, if the EDM is not used above
                if (~Zmn_calculated)
                    num_std_zmn_entries_calculated = num_std_zmn_entries_calculated + 1;

                    %fprintf('Processing element %d,%d\n',mm,nn);
                    % There are four terms to add here;
                    % From integrals over source faces n+ & n- evaluated at centres of field faces m+ and m-.
                    % The n integrals associate with q, the m, with p.
                    % In datastructures EDGECONXELEMS & DOFLOCALNUM, second index 1
                    % associates with +, 2 with -. 

                    %qq_pls = EDGECONXELEMS(nn,1);
                    qq_pls = Solver_setup.rwg_basis_functions_trianglePlus(nn);
                    %qq_mns = EDGECONXELEMS(nn,2);
                    qq_mns = Solver_setup.rwg_basis_functions_triangleMinus(nn);
                
                    triangle_tn_plus_free_vertex = Solver_setup.rwg_basis_functions_trianglePlusFreeVertex(nn);
                    triangle_tn_minus_free_vertex = Solver_setup.rwg_basis_functions_triangleMinusFreeVertex(nn);
                    
                    % First, find contribution from n+ and n- faces evaluated at m+
                                    
                    % --------------------------------------------------------------                
                    % Look at Tm+
                    % --------------------------------------------------------------
                    % -- Contribution from Tn+
                    
                    [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_pls,qq_pls,mm,nn,triangle_tn_plus_free_vertex,...
                        k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                    
                    % [DBD2011] implementation below. I think there is a sign issue here.
                    %Amn_pls_source_pls = MagVecPot;
                    %Phi_mn_pls_source_pls = -ScalPot;
                    % [DL - 2018] implementation: Swop the signs around
                    Amn_pls_source_pls = -MagVecPot;
                    Phi_mn_pls_source_pls = +ScalPot;

                    
                    % -- Contribution from Tn-
                    
                    [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_pls,qq_mns,mm,nn,triangle_tn_minus_free_vertex,...
                        k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                    % [DBD2011] implementation below. I think there is a sign issue here.
                    %Amn_pls_source_mns = - MagVecPot;
                    %Phi_mn_pls_source_mns = +ScalPot;
                    % [DL - 2018] implementation: Swop the signs around
                    Amn_pls_source_mns = + MagVecPot;
                    Phi_mn_pls_source_mns = -ScalPot;
                                    
                    Amn_pls = Amn_pls_source_pls + Amn_pls_source_mns;
                    Phi_mn_pls = Phi_mn_pls_source_pls + Phi_mn_pls_source_mns;
                    
                    % --------------------------------------------------------------                
                    % Look at Tm-
                    % --------------------------------------------------------------
                    % -- Contribution from Tn+
                    [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_mns,qq_pls,mm,nn,triangle_tn_plus_free_vertex,...
                        k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                    % [DBD2011] implementation below. I think there is a sign issue here.
                    %Amn_mns_source_pls = MagVecPot;
                    %Phi_mn_mns_source_pls = -ScalPot;
                    % [DL - 2018] implementation: Swop the signs around
                    Amn_mns_source_pls = -MagVecPot;
                    Phi_mn_mns_source_pls = +ScalPot;
                    
                    % -- Contribution from Tn-                
                    [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell,pp_mns,qq_mns,mm,nn,triangle_tn_minus_free_vertex,...
                        k,r_c,quad_pts,sing,eps_0,mu_0,omega);
                    % [DBD2011] implementation below. I think there is a sign issue here.
                    %Amn_mns_source_mns = - MagVecPot;
                    %Phi_mn_mns_source_mns = +ScalPot;
                    % [DL - 2018] implementation: Swop the signs around
                    Amn_mns_source_mns = +MagVecPot;
                    Phi_mn_mns_source_mns = -ScalPot;
                    
                    Amn_mns = Amn_mns_source_pls + Amn_mns_source_mns;
                    Phi_mn_mns = Phi_mn_mns_source_pls + Phi_mn_mns_source_mns;
                           
                    % Assemble with eq. 17 in [RWG82]
                    Z.values(mm,nn) = 1i*omega*...
                        (dot(Amn_pls',rho_c_pls(mm,:))/2 + dot(Amn_mns',rho_c_mns(mm,:))/2) + Phi_mn_mns - Phi_mn_pls;
                    
                    %mm
                    %nn
                    Z.values(mm,nn) = ell(mm)* Z.values(mm,nn);                    
                end % if (~Zmn_calculated)
                
            end % for nn = 1:NUM_DOFS
        end %for mm = 1:NUM_DOFS
    end %for freq_index = 1:Solver_setup.frequencies.freq_num    
    
    if (Const.useEDM)
        EDM_mat_perc = num_edm_approximations / (num_dofs*num_dofs) * 100;
        message_fc(Const,sprintf('  Number of EDM approximations   : %d (%f perc. of MoM matrix)',num_edm_approximations, EDM_mat_perc));
        message_fc(Const,sprintf('  Number of Standard Zmn entries : %d ',num_std_zmn_entries_calculated));
    end%if
    
end %function FillZMatrixByEdge

% =================================================================================
% Local function for evaluating the potential integrals, as done in [1]
% =================================================================================
function [MagVecPot,ScalPot] = Potentials(elements,node_coord,ell, field_pt,source_pt,field_edge, source_edge, ii, k,r_c,quad_pts,sing,eps_0,mu_0,omega)
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

    % [RWG82, Eq. (32) - without sign
    MagVecPot = mu_0*ell(source_edge)/(4*pi)*...
        ( r(1,:)*Ipq_xi + r(2,:)*Ipq_eta + r(3,:)*Ipq_zeta - rii(1,:)*Ipq);

    % [RWG82, Eq. (33) - without sign
    ScalPot = ell(source_edge)/(1i*2*pi*omega*eps_0) * Ipq;
        
end % function Potentials

