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

    % -- Calculate the 
    [zMatrices] = FillZMatrixByEdge(Const, Solver_setup);

    % -- Calculate the yVectors, i.e. the RHS
    yVectors = []; % TO-DO: calculate
    