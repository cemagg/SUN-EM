function [xySamples] = extractDGFMcircularInterpolationGrid(Const, numArrayEls, xyPoints)
    %extractDGFMcircularInterpolationGrid
    %   Usage:
    %       [xySamples] = extractDGFMcircularInterpolationGrid(Const, xyPoints)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       xyPoints
    %           A 2D matrix containing (as rows) the X and Y positions of each array element in the array  
    %
    %   Output Arguments:
    %       xySamples
    %           Vector containing the sampling indices for each of the
    %           array elements, i.e. where the full/reference answer will
    %           be calculated - in this case, e.g. the i-DGFM solution.
    %
    %   Description:
    %       Given a set of potentially random array element positions (XY)
    %       laid out in a circular grid, we return a set of more or less
    %       equally spaced indices representing a uniformly sampled subset
    %       of these elements. Note: This currently assumes a 2D circular 
    %       irregular array layout.

    narginchk(3,3);
    
    % Local debug flag
    LOCAL_DEBUG = true;
    
    % Set up a random sampling grid now on the array lattice
    numberOfsamples = ceil(Const.DGFMinterpolationSamplingFactor * numArrayEls);
    
    % TO-DO: Here we would need to do a bit of additional work to get the
    % sampling grid better set-up (i.e. more uniformly). Look at e.g. a 
    % Latin Hypercube - or, just randomly sample the 4 quadrants
    % separately.
    xySamples = randsample(numArrayEls, numberOfsamples);
    
    % Plot the array + the DGFM sampling positions (should be equally distributed)
    if (LOCAL_DEBUG)
        figure;
        hold on;
        grid on;
        plot(xyPoints(:,1),xyPoints(:,2), 'o', 'MarkerSize',10, 'LineWidth', 2);
        plot(xyPoints(xySamples(:),1),xyPoints(xySamples(:),2), 'o', 'MarkerFaceColor','r','MarkerSize',10);
        %legend('Array layout', 'Sampling point', 'FontSize', 12);
        xlabel('Distance [m]', 'FontSize',12);
        ylabel('Distance [m]', 'FontSize',12);
        
        % display also the element numbers 
        for arrayIndex = 1:numArrayEls
            text(xyPoints(arrayIndex,1), xyPoints(arrayIndex,2), num2str(arrayIndex),'FontSize',16);
        end%for
    end%if