function [Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const)
    %extractFEKOMoMmatrixEq v0.1
    %   Date: 28.11.2013
    %   Usage:
    %           [Const, zMatrices, yVectors, xVectors] = extractFEKOMoMmatrixEq(Const)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containging settings of which solver to run,
    %           as well as initial basis function setup
    %
    %   Output Arguments:
    %       Const
    %           Struct containing basis function setup
    %       zMatrices
    %           The Z-matrices data
    %       yVectors
    %           The Yrhs-vector data
    %       xVectors
    %           The Xsol-vector data (i.e. MoM solution of FEKO)
    %
    %   Description:
    %       Extracts the FEKO MoM matrix from the FEKO *.out, *.mat, *.str and *.rhs files
    %
    %   =======================
    %   Written by Danie Ludick on November 28, 2013.
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    narginchk(1,1);

    % -- Read Z
    [zMatrices] = readFEKOZMatrixFromFile(Const, Const.FEKOmatfilename);

    % -- Read X (FEKO MoM solution)
    [xVectors] = readFEKOXvectorFromFile(Const, Const.FEKOstrfilename);

    % -- Read Y (first check that the dimensions of Zmn and Xsol correspond)
    if ( (zMatrices.mBasis == zMatrices.nBasis) && (zMatrices.mBasis == xVectors.numMoMbasis) )
        Const.numMoMbasis = zMatrices.mBasis; % Nrwg
    else
        message_fc(Const,'Invalid matrix equation dimensions');
        error ('Invalid matrix equation dimensions');
    end
    [yVectors] = readFEKOYvectorFromFile(Const, Const.FEKOrhsfilename);
    
    % Consistency check - the number of RHS vectors in yVectors and the number of reference
    % solutions in xVectors should be the same
    if (xVectors.numSols ~= yVectors.numRhs)
        message_fc(Const,sprintf('xVectors.numSols: %d, yVectors.numRhs: %d', xVectors.numSols, yVectors.numRhs));
        message_fc(Const,'The number of solutions does not agree with the number of excitation vectors');
        error ('Invalid matrix equation dimensions');
    end

    % Initialisation of FEKO-Connect based on input parameter settings in Const.
    Const = sunem_init(Const, yVectors);
