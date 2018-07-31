function [Const, zMatrices, yVectors] = extractSUNEMMoMmatrixEq(Const, Solver_setup)
    %extractSUNEMMoMmatrixEq
    %   Date: 2018.06.10
    %   Usage:
    %           [Const, zMatrices, yVectors] = extractSUNEMMoMmatrixEq(Const, Solver_setup)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containging settings of which solver to run,
    %           as well as initial basis function setup
    %       Solver_setup
    %           The solution setup, i.e. geometry, basis function setup, etc.
    %
    %   Output Arguments:
    %       Const
    %           Struct containing basis function setup
    %       zMatrices
    %           The Z-matrices data calculated internally
    %       yVectors
    %           The Yrhs-vector data calculated internally
    %
    %   Description:
    %       Extracts the FEKO MoM matrix from the FEKO *.out, *.mat, *.str and *.rhs files
    %
    %   =======================
    %   Written by Danie Ludick on 2018.06.10
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za
    %   Credit to: Prof. David B. Davidson for the FillZMatrixByEdge routine that is based on his
    %   MATLAB implementation as detailed in [1]
    %   
    %   References: David B. Davidson, Computational Electromagnetics for RF and Microwave Engineering, 
    %               Second Edition, (see Chapter 6)

    narginchk(2,2);

    message_fc(Const,' ');
    message_fc(Const,...
        '------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Calculating MoM matrix internally'));
    if (Const.useEDM)
        message_fc(Const,sprintf('  EDM active'));
    end%if
    if (Const.use_CPP_engine)
        message_fc(Const,sprintf('  Using C++ engine'));
    end%if

    % -- Calculate the Z-matrix (edge-based fill)
    zMatfilltime=tic;
        if (Const.use_CPP_engine)
            %TO-DO: Add here the fast C++ Z matrix engine (faster)
            %[zMatrices] = FillZMatrixByEdgeCPP(Const, Solver_setup);
            message_fc(Const,'  C++ engine not yet active');
            error('  C++ engine not yet active');
        else
            [zMatrices] = FillZMatrixByEdge(Const, Solver_setup);
        end%if
    totMatrixSetupTime = toc(zMatfilltime);

    % Define a normally incident, X-directed plane wave:
    EMag = 1;
    theta_0 = 0;
    phi_0 = 0;
    
    % -- Calculate the yVectors, i.e. the RHS
    vVecfilltime=tic;        
        yVectors.values = FillVVector(Const,Solver_setup,EMag,theta_0,phi_0);        
    totRHSvecSetupTime = toc(vVecfilltime);
    
    % Output total time
    message_fc(Const,sprintf('Time for Z-matrix setup       : %f sec.',totMatrixSetupTime));
    message_fc(Const,sprintf('Time for V-vector (RHS) setup : %f sec.',totRHSvecSetupTime));
    