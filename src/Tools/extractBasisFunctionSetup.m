function Const = extractBasisFunctionSetup(Const, yVectors)
    %extractBasisFunctionSetup
    %   Usage:
    %           extractBasisFunctionSetup(Const)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containging settings of which solver to run,
    %           as well as initial basis function setup
    %       yVectors
    %           The Yrhs-vector data
    %
    %   Output Arguments:
    %       Const
    %           Struct containing basis function setup
    %
    %   Description:
    %       Extracts basis function setup and array numbering, etc. as was
    %       read / parsed from the FEKO *.out, *.mat, *.str and *.rhs files
    %
    %   =======================
    %   Written by Danie Ludick on November 24, 2013.
    %   Last updated on November 28, 2013.
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    error(nargchk(2,2,nargin));

    % LEGACY CODE BELOW --> WE NEED TO REFACTOR THIS

    % Initialise the total BFs associated with the finite array geometry
    % Note, if the NGF is also used, then we need to discard the static domain
    %Const.numArraybasis = Const.numMoMbasis - Const.numNGFbasis;
    Const.numArraybasis = Const.numMoMbasis;
    Const.numMoMbasisPerElement = Const.numArraybasis/Const.numArrayElements;