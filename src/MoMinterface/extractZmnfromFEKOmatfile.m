function [Zmn] = extractZmnfromFEKOmatfile(Const, zMatrices, freq, observBFs, sourceBFs)
    %buildMoMblock
    %   Usage:
    %       [Zmn] = extractZmnfromFEKOmatfile(Const, zMatrices, observBFs, sourceBFs)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       zMatrices
    %           The Z-matrices data (as read by FEKO from the *.mat file)
    %       freq
    %           The frequency index for which this MoM matrix should be extracted
    %       observBFs
    %           The observer basis function indices (a vector)
    %       sourceBFs
    %           The source basis function indices (a vector)
    %
    %   Output Arguments:
    %       Zmn
    %           A submatrix corresponding to the observer and source BFs,
    %           as extracted from zMatrices (FEKO *.mat file)
    %
    %   Description:
    %       Extracts a specific submatrix from zMatrices, i.e. as read from the FEKO
    %       *.mat file (use a delay here to "simulate" a slower Amn calculation). 
    %       This is done by the flag Const.fastBuildMoMblock.
    %
    %   =======================
    %   Written by Danie Ludick on March 24, 2014.
    %   Last updated on June 25, 2017.
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    error(nargchk(5,5,nargin));

    if (Const.fastBuildMoMblock)
        % Fast algorithm -- Extract the matrix elements using a vectorised operation
        Zmn = zMatrices.values(observBFs, sourceBFs, freq);
    else
        % Slow algorithm -- Extract the matrix elements one-by-one (causes a delay)
        Zmn = complex(zeros(length(observBFs),length(sourceBFs)));
        % We need to implement a loop here (i.e. element per element)
        for m=1:length(observBFs)
            for n=1:length(sourceBFs)
                Zmn(m,n) = zMatrices.values(observBFs(m), sourceBFs(n), freq);
                % Removed the following as the delay resulted in too long
                % runtimes (Not needed as the element-by-element access here is
                % slow enough):
                % Pause for 1ms to simulate calculating Amn
                % pause(1e-20);
            end%for
        end%for
    end%if
