function [Zmn,U,V] = calcZmn(Const, zMatrices, freq, m, n, observBFs, sourceBFs)
    %calcZmn
    %   Usage:
    %       [zMatrices] = calcZmn(Const, zMatrices, m, n, observBFs, sourceBFs)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       zMatrices
    %           The Z-matrices data
    %       freq
    %           The frequency index for which this MoM matrix should be extracted
    %       m,n
    %           The matrix row and column subscripts
    %
    %   Output Arguments:
    %       Zmn
    %           The coupling matrix Zmn that was calculated with the MoM (fast or slow algorithm).
    %           TO-DO: Future extension here can perhaps construct the MoM matrix ourselves.
    %
    %   Description:
    %       Calculates the coupling matrix Zmn.
    %
    %   =======================
    %   Written by Danie Ludick on March 24, 2014.
    %   Last updated on June 25, 2017.
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    error(nargchk(7,7,nargin));

    Zmn = extractZmnfromFEKOmatfile(Const, zMatrices, freq, observBFs, sourceBFs);

    % TO-DO: Add here more algorithms as we develop them.
    U = [];
    V = [];
