function [V] = FillVVector(Const, Solver_setup, EMag,theta_0,phi_0)
    %FillVVector
    %   Date: 2018.06.12
    %   Usage:
    %       [V] = FillVVector(EMag,theta_0,phi_0)
    %
    %   Input Arguments:
    %       Const: 
    %           A global struct containing:
    %       Solver_setup
    %           A global struct containing solver setup details, e.g. frequency range,
    %           geometry data, basis function setup, etc.
    %   Output Arguments:
    %       V
    %           The Y-vector (RHS) for a normally incident, X directed
    %           plane wave
    %
    %   Description:
    %       Fills the Y-vector, i.e. RHS in the MoM matrix equation for a
    %       normally incident, X-directed plane wave. 
    %       TO-DO: See eq 22 and 23 [RWG82] for other PW incident angle
    %       definitions.
    %
    %   =======================
    %   Written by Danie Ludick on 2018.06.12
    %   Stellenbosch University
    %   Email: dludick.sun.ac.za 
    %   Credit to: Prof. David B. Davidson for the FillVVector routine that is based on his
    %   MATLAB implementation as detailed in [1]
    %   
    %   References: 
    %   [1] David B. Davidson, Computational Electromagnetics for RF and Microwave Engineering, 
    %       Second Edition, (see Chapter 6)

    message_fc(Const, sprintf('  Calculating V-vector (RHS) '));
    
    num_dofs = Solver_setup.num_metallic_edges;       % Replacing global NUM_DOFS
    ell = Solver_setup.rwg_basis_functions_length_m;  % Replacing global ELL
    rho_c_pls = Solver_setup.rho_c_pls;
    rho_c_mns = Solver_setup.rho_c_mns;

    V = zeros(num_dofs,1); % RHS matrix
    E = zeros(num_dofs,3); % Sampled incident electric field
    % x-polarized field, unity magnitude
    E(:,1) = EMag; % Special case
    %E(:,2) = 1; % Special case

    % We only allow normally incident, x-directed plane wave at the moment
    if (theta_0 ~= 0 || phi_0 ~= 0 )
        message_fc(sprintf('Only normally incident plane wave implemented at present '));
        error('Only normally incident plane wave implemented at present ')
    end

    for mm = 1:num_dofs
        %pp_pls = EDGECONXELEMS(mm,1);            
        pp_pls = Solver_setup.rwg_basis_functions_trianglePlus(mm);
        %pp_mns = EDGECONXELEMS(mm,2);
        pp_mns = Solver_setup.rwg_basis_functions_triangleMinus(mm);
        if (mm == 894)
            a = 1;
        end%if
        % TO-DO: Danie, check the following, I uncommented now the
        % E(pp_pls,:) access, as it resulted in a segmentation violation. I
        % think we are sampling it on the edge, mm.
        %V(mm) = dot(E(pp_pls,:)',rho_c_pls(mm,:))/2 + dot(E(pp_mns,:)',rho_c_mns(mm,:))/2;
        % or, the incident electric field is not correctly defined (i.e.
        % sized above).
        V(mm) = dot(E(mm,:)',rho_c_pls(mm,:))/2 + dot(E(mm,:)',rho_c_mns(mm,:))/2;
        V(mm) = ell(mm)*V(mm);
    end


