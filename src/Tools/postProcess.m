function [results] = postProcess(Const, solution)
    %postProcess
    %   Usage:
    %       [results] = postProcess(Const, results
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containging settings of which solver to run
    %       solution
    %           The solution struct, e.g. that pertaining to the MoM/HARP, etc.
    %           solution
    %
    %   Output Arguments:
    %      solution.
    %           Structs containing results, e.g. the far-field pattern
    %
    %   Description:
    %       Runs different post-processing utilities, e.g. calculate the directivity
    %       pattern.

    error(nargchk(2,2,nargin));

    % -- Directivity pattern
    if (Const.calcDirectivity)
        results = calculateDirectivity(Const, solution);
    end%if
