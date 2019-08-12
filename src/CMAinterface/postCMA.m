function [cma] = postCMA (Const, cma, FEKO_mode_data)
    %postCMA v0.1
    %   Date: 26.05.2014
    %   Usage:
    %       postCMA (Const, cma)
    %
    %   Input Arguments:
    %       Const   
    %           A global struct containg the parameters, e.g., debug
    %           settings, *.str filename, etc.
    %       cma
    %           Struct containing the CMA results (e.g. the eigen values and 
    %           eigenvectors) as calculated in runCMAsolver.m
    %       FEKO_mode_data
    %           FEKO's CMA data that is extracted from the *.out file
    %
    %   Output Arguments:
    %       
    %   Description:
    %       Plots various results (e.g. the eigenvalues vs. frequency using
    %       both tracked and untracked mode ranking values), the mode
    %       orthogonality, etc. See also use of Const.postCMAplots (works
    %       with binary format)
    %
    %       TO-DO: Enable the following:
    %       Const.postCMAplots  = 1: Eigenvalue vs. Frequency
    %                             2: Char. angle vs. Frequency
    %                             4: Mode Orthogonality vs. Frequency
    %                             Combinations, e.g.
    %                             3: Do (1) and (2)
    %                             7: All of the above
    %
    %   =======================
    %   Written by Danie Ludick on May 26, 2014
    %   Last updated on May 26, 2014.
    %   EM Software & Systems - S.A. (Pty) Ltd.
    %   Email: dludick@emss.co.za
    
    %   References:
    %   ----------
    
    error(nargchk(3,3,nargin));
    
    % Close all previous images
    %close all;
    
    message_fc(Const,sprintf(''));

    % Initialise variables
    numModes = Const.numCMAmodes;
    numFreq =  cma.numFreq;
    
    % Note: When including now the CMA for finite array analysis, i.e. 
    % by applying also the active impedance matrix (Zact), we will be
    % analysing each array element individually. The index of the element
    % is stored in 
    
    message_fc(Const,sprintf('Post processing CMA for array element: %d',Const.arrayID));
    
    % ===============================
    % TO-DO: This should also be set from input file (see above comment)
    plotEigvalVsFreq    = false;
    plotCharAngleVsFreq = true;
    if (Const.checkRXorth)
        plotModeOrthVsFreq = true;
    else
        plotModeOrthVsFreq = false;
    end%if        
    
    plotXsolVsFreq = false;       % Plots the eigenvectfor coefficients vs. frequency
    
    % -------------------------------
    % For the adaptive extrapolation:
    % -------------------------------
    plotExtrSamplesVsFreq.val = false; % Plots the extrpolated eigenvalue / charact. angle
                                      % curve based on the sorted values superimposed ON the
                                      % unsorted corresponding curve

    plotExtrSamplesVsFreq.type = 0;   % 0 : Eigen Val.
                                      % 1 : Charact. Angle.
    
    plotExtrSamplesVsFreq.mode = 4;            % The mode index used for the extrapolation

    plotExtrSamplesVsFreq.numFreqSamples = 2; % The number of frequency samples to use (we will then
                                               % extrapolate to this number + 1)

    plotExtrSamplesVsFreq.extrMethod = 1;      % 0: Robert's method (in ADAPTFEKO)
                                               % 1: MATLAB interp1 routine
                                               %   (see also different
                                               %   subsettings there, e.g.
                                               %   'linear', 'spline', etc.
    % ------------
                                      
    % To avoid confusion when plotting the extrapolated sampling values, 
    % enable only the Eig. Vals vs. Frequency (untracked) / Charact. Angle vs. Freq 
    % so that we can have an idea of where the extrapolation lies
    if (plotExtrSamplesVsFreq.val)
        if (plotExtrSamplesVsFreq.type == 0)
            % -- Eigenvalue vs. Freq. curve
            plotEigvalVsFreq    = true;            
            plotCharAngleVsFreq = false;
        else
            % -- Charact. Angle vs. Freq. curve
            plotEigvalVsFreq    = false;
            plotCharAngleVsFreq = true;
        end
        
        % -- The following is false regardless
        plotModeOrthVsFreq = false;
        plotXsolVsFreq = false;
    end%if
                                         
    modeNum = 15; % Plot the first modeNum modes
    
    % Use FEKO's modal ranks for testing
    useFEKOmodrank = false;

    % Check to see whether we used the FEKO data. If not, then the above is
    % set to false regardless:
    if (isempty(FEKO_mode_data))
        useFEKOdata = false;
        useFEKOmodrank = false;
    end%if
    
    % Calculate the correct frequency range for the plots
    % Already done in the driver - just do an internal consistency
    % check here
    freqData = Const.freqData;
    if (length(freqData)~=numFreq)
        message_fc(Const,sprintf('[postCMA] Internal Error 1.0'));
        error(['[postCMA] Internal Error 1.0']);
    end%if
    
    if (useFEKOdata)
        % ===============================
        % Compare now FEKO's eigenvalues with that of FEKCMA (all frequencies,
        % all requested modes)  
        for freq = 1:numFreq
            % NOTE: The following loop only looks over the number of requested
            % eigenvalues, as this is what could be read from the *.out file
            for mm = 1:FEKO_mode_data.num_requested_modes
                err = calculateErrorNormPercentage( abs(FEKO_mode_data.eigr(mm,freq)+1i*FEKO_mode_data.eigi(mm,freq)), ...
                                                    abs(cma.EigVal(mm,freq)) );
                % Check whether these values are more or less the same
                if (err >= 1e-1)
                    message_fc(Const,sprintf('[postCMA] Eigen val. %d at Freq. inc. %d differs too much from FEKO',mm,freq));
                    message_fc(Const,sprintf('[postCMA] FEKCMA mode-val: %.5E',cma.EigVal(mm,freq)));
                    message_fc(Const,sprintf('[postCMA] FEKO mode-val: %.5E + j*%.5E',FEKO_mode_data.eigr(mm,freq),FEKO_mode_data.eigi(mm,freq)));
                    % Only error when we do not explicitely symmetrise the Z
                    % matrix here in FEKO (otherwise we will get slightly
                    % different answers):
                    if (~Const.symZmat)
                        error(['[postCMA] Internal Error 2.0']);
                    end%if
                end%if
            end
        end
    end%if (useFEKOdata)
    
    % ===============================
    % Extract some mode tracking data here
    if (Const.trackModes)
        
        % Remove the following -- irretating if we do debugging and limit
        % the number of modes reported by FEKO:
        if (false && (useFEKOdata))
            % Do a comparison between FEKO and FEKCMA's ranking values:
            cma.modetrack_alg1_diff = FEKO_mode_data.mode_ranks - cma.mode_ranks;
            if (~isempty(find(FEKO_mode_data.mode_ranks - cma.mode_ranks)))
                message_fc(Const,sprintf('[postCMA] Difference between FEKO and FEKCMA tracking alg. 1'));
            end
            
            % Check adaptive ranking values
            cma.modetrack_alg2_diff = FEKO_mode_data.mode_ranks - cma.mode_ranks_adaptive;
            if (~isempty(find(FEKO_mode_data.mode_ranks - cma.mode_ranks_adaptive)))
                message_fc(Const,sprintf('[postCMA] Difference between FEKO and FEKCMA tracking alg. 2'));
            end
        end
        
        % Set here the mode ranks to use when plotting the tracked results.
        if (Const.trackModes == 1)
            % -- Standard tracking
            modrank = cma.mode_ranks;
        elseif (Const.trackModes == 2)
            % -- Adaptive tracking
            modrank = cma.mode_ranks_adaptive;    
        end%if
        
    end%if
    
    % ============================================================
    %              Untracked Eigenvalue vs. Frequency
    % ============================================================
    if (plotEigvalVsFreq)
        figure
        hold on;
        grid on;
        box on;
        ha = gca; % get current axes handle
        set(ha, 'FontSize', 16);
        set(ha, 'LineWidth', 2);
        xlabel('Frequency [GHz]');
        ylabel('Eigen Value');
        title('Untracked Eigenvalue vs. Frequency');
        for mm = 1:modeNum
            col = rand(1,3);
            plot(freqData,cma.EigVal(mm,:),'-o','color',col,'LineWidth',2.5,...
                'MarkerSize',9,'MarkerFaceColor',col);
            legendInfo{mm} = ['Mode ' num2str(mm)]; 
        end%for
        legend(legendInfo);
    end%if (plotEigvalVsFreq)
    
    % ============================================================
    %              Char. Angle vs. Frequency
    % ============================================================
    
    if (plotCharAngleVsFreq)
        
        % ------------------------------------
        %              Untracked
        % ------------------------------------
        
        figure
        hold on;
        grid on;
        box on;
        ha = gca; % get current axes handle
        set(ha, 'FontSize', 16);
        set(ha, 'LineWidth', 2);
        xlabel('Frequency [GHz]');
        ylabel('Characteristic Angle');
        
        if (Const.useZactCMA)
            title(sprintf('Untracked Characteristic Angle vs. Frequency for Array El. %d', ...
                Const.arrayID));
        else
            title('Untracked Characteristic Angle vs. Frequency');
        end
        
        % Plot here all the untracked modes (indices limited to number
        % calculated internally)        
        for mm = 1:Const.numModesCalc % Originally Const.numModesCalc
            % Remove the following - plot all indices
            if (true)
                if (mm > modeNum)
                    % Plot only the number of modes requested
                    break;
                end%if
            end%if
            col = rand(1,3);
            charAngle = 180 - atan(cma.EigVal(mm,:))*180/pi;
            plot(freqData,real(charAngle),'-o','color',col,'LineWidth',2.5,...
                'MarkerSize',9,'MarkerFaceColor',col);
            legendInfo{mm} = ['Mode ' num2str(mm)];
        end%for
        legend(legendInfo);
                
        if (~plotExtrSamplesVsFreq.val && Const.trackModes)
            % For extrapolation - plot only the above untracked curve
            
            if (useFEKOmodrank && useFEKOdata)
                % ------------------------------------
                %              Tracked (FEKO results)
                % ------------------------------------
                % Use FEKO's mode ranking vector
                
                figure
                hold on;
                grid on;
                box on;
                ha = gca; % get current axes handle
                set(ha, 'FontSize', 16);
                set(ha, 'LineWidth', 2);
                xlabel('Frequency [GHz]');
                ylabel('Characteristic Angle');
                title('Tracked (FEKO) Characteristic Angle vs. Frequency');
                
                % First calculate the total number of curves that has been logged
                % The number of which might be more than the # of modes requested, 
                % as some modes might stop while others (new indices) start
                maxNumTraces = max(max(FEKO_mode_data.mode_ranks));
                message_fc(Const,sprintf('[postCMA] Maximum number of traces available: %d',maxNumTraces));
                % For each of the trace ID's extract the correct data (eigen value)
                % at each of the frequencies
                for traceID = 1:maxNumTraces
                    
                    % Remove the following - plot all indices
                    if (true)
                        if (traceID > modeNum)
                            % Plot only the number of modes requested
                            break;
                        end%if
                    end%if
                    
                    % Loop over each of the frequencies and extract the correct trace                    
                    % info (eigenvalue and frequency values)
                    trackedEigVal = [];
                    traceFreqData = [];
                    for freqID = 1:numFreq
                        % Extract the eigenvalue data associated with traceID
                        pos = find(FEKO_mode_data.mode_ranks(:,freqID)==traceID);
                        if (~isempty(pos))
                            traceFreqData = [traceFreqData freqData(freqID)];
                            trackedEigVal = [trackedEigVal cma.EigVal(pos,freqID)];
                        end
                    end%for
                    
                    %Plot the resulting (tracked) characteristic angle here
                    col = rand(1,3);
                    charAngle = 180 - atan(trackedEigVal)*180/pi;
                    % Take note here to use the trace's frequency range
                    plot(traceFreqData,charAngle,'-o','color',col,'LineWidth',2.5,...
                        'MarkerSize',9,'MarkerFaceColor',col);
                    legendInfo{traceID} = ['Mode ' num2str(traceID)];            
                end%for        
                legend(legendInfo);
            end%if
            
            % ------------------------------------
            %              Tracked
            % ------------------------------------
            
            figure
            hold on;
            grid on;
            box on;
            ha = gca; % get current axes handle
            set(ha, 'FontSize', 16);
            set(ha, 'LineWidth', 2);
            xlabel('Frequency [GHz]');
            ylabel('Characteristic Angle');
            
            if (Const.useZactCMA)
                title(sprintf('Tracked Characteristic Angle vs. Frequency for Array El. %d', ...
                Const.arrayID));
            else
                title('Tracked Characteristic Angle vs. Frequency');
            end
                        
            % First calculate the total number of curves that has been logged
            % The number of which might be more than the # of modes requested, 
            % as some modes might stop while others (new indices) start
            maxNumTraces = max(max(modrank));
            message_fc(Const,sprintf('[postCMA] Maximum number of traces available: %d',maxNumTraces));
            % For each of the trace ID's extract the correct data (eigen value)
            % at each of the frequencies
            for traceID = 1:maxNumTraces        
               % Remove the following - plot all indices
                if (true)
                    if (traceID > modeNum)
                        % Plot only the number of modes requested
                        break;
                    end%if
                end%if
                % Loop over each of the frequencies and extract the correct trace
                % info (eigenvalue and frequency values)
                trackedEigVal = [];
                traceFreqData = [];
                for freqID = 1:numFreq
                    % Extract the eigenvalue data associated with traceID
                    pos = find(modrank(:,freqID)==traceID);
                    if (~isempty(pos))
                        traceFreqData = [traceFreqData freqData(freqID)];
                        % DJdbg --> remove
                        if ((mm == 5) && (freqID == 11))
                            freqID
                            a = 1;
                        end%if
                        trackedEigVal = [trackedEigVal cma.EigVal(pos,freqID)];
                    end
                end%for
                
                %Plot the resulting (tracked) characteristic angle here
                col = rand(1,3);
                charAngle = 180 - atan(trackedEigVal)*180/pi;
                % Take note here to use the trace's frequency range
                plot(traceFreqData,real(charAngle),'-o','color',col,'LineWidth',2.5,...
                    'MarkerSize',9,'MarkerFaceColor',col);
                legendInfo{traceID} = ['Mode ' num2str(traceID)];            
            end%for        
            legend(legendInfo);
            
        end%if (~plotExtrSamplesVsFreq.val)  
    end%if (plotCharAngleVsFreq)

    % ============================================================
    %              Mode Orthogonality, etc.
    % ============================================================

    if (plotModeOrthVsFreq)
        figure
        hold on;
        grid on;
        box on;
        ha = gca; % get current axes handle
        set(ha, 'FontSize', 16);
        set(ha, 'LineWidth', 2);
        xlabel('Frequency [GHz]');
        ylabel('Number of R and X orthogonal modes');
        title('R and X orthogonality of modes vs. Frequency');
        % FEKCMA data
        plot(freqData,cma.numOrthModes,'-o','color','k','LineWidth',2.5,...
            'MarkerSize',9,'MarkerFaceColor','k');
        if (useFEKOdata)
            % FEKO data
            plot(freqData,FEKO_mode_data.num_orthogonal_modes,'-o','color','b','LineWidth',2.5,...
                'MarkerSize',9,'MarkerFaceColor','b');
        end%if
        if (Const.symZmat)
            fekcmastr = 'FEKCMA (symm. [Z])';
        else
            fekcmastr = 'FEKCMA';
        end%if
        
        if (useFEKOdata)
            legend(fekcmastr,'FEKO');
        else
            legend(fekcmastr);
        end%if
    end%if (plotModeOrthVsFreq)
    
    % For the tracking we would like to use ADAPTIVE interpolation to 
    if (plotXsolVsFreq)
        
        figure
        hold on;
        grid on;
        box on;
        ha = gca; % get current axes handle
        set(ha, 'FontSize', 16);
        set(ha, 'LineWidth', 2);
        xlabel('Frequency [GHz]');
        ylabel('Xi');
        
        XsolModeIndx = 1;
        title(sprintf('Mode %d expansion coefficients vs. Frequency',XsolModeIndx));
        
        % Sample here the (XsolSamples)-expansion coefficients of mode 1
        XsolSamples = 5;
        legendInfo = [];
        for xi = 1:XsolSamples
            col = rand(1,3);
            % Extract the calculated eigenvector coefficients
            xsol = zeros(numFreq,1);
            for freqID = 1:numFreq
                xsol(freqID,1) = abs(cma.Isol(XsolModeIndx,xi,freqID));
            end%for                
            plot(freqData,xsol,'-o','color',col,'LineWidth',2.5,...
                    'MarkerSize',9,'MarkerFaceColor',col);
            legendInfo{xi} = ['Xi :' num2str(xi)];
        end%for
        legend(legendInfo);
        
    end%if
    
    % ------------------------------------
    % 29.05.2014: What would happen if we extrapolate the eigenvalue /
    % charact. angle curves to get our match? Plot here the extrapolation
    % of a single (sorted) curve on the unsorted plots.
   
    if (plotExtrSamplesVsFreq.val && Const.trackModes)
        
        %figure; % Results added to untracked eigenvalue curve
        hold on;
        grid on;
        box on;
        ha = gca; % get current axes handle
        set(ha, 'FontSize', 16);
        set(ha, 'LineWidth', 2);
        xlabel('Frequency [GHz]');
        
        modeIndx = plotExtrSamplesVsFreq.mode; % Mode to be extrapolated
        if (plotExtrSamplesVsFreq.type == 0)
            % -- Eigen Value
            title(sprintf('Extrapolated Eigen Val. of mode %d vs. Frequency',modeIndx));
        else
            % -- Characteristic Angle
            title(sprintf('Extrapolated Charact. Angle of mode %d vs. Frequency',modeIndx));
        end%if
        
        col = rand(1,3);
        
        trackedVal = [];
        traceFreqData = [];
        for fID = 1:plotExtrSamplesVsFreq.numFreqSamples
            % Extract the tracked data associated with modeIndx
            pos = find(cma.mode_ranks_adaptive(:,fID)==modeIndx);
            if (~isempty(pos))
                traceFreqData = [traceFreqData Const.freqData(fID)];
                if (plotExtrSamplesVsFreq.type == 0)
                    % -- Eigen Value
                    trackedVal = [trackedVal cma.EigVal(pos,fID)];
                else
                    % -- Characteristic Angle
                    trackedVal = [trackedVal cma.CharAngle(pos,fID)];                        
                end%if
            else
                % Error for now - extrapolate select curves that will work
                message_fc(Const,sprintf('[postCMA] Internal Error 3.0'));
                error(['[postCMA] Internal Error 3.0']);
            end%if
        end%for
        
        % Create the frequency samples at which we would like our extrapolation
        freqExtr = [traceFreqData Const.freqData(plotExtrSamplesVsFreq.numFreqSamples+1)].';
        
        if (plotExtrSamplesVsFreq.extrMethod==0)
            % -- Use Robert's adaptfeko interpolation routines            
            [ValExtr] = extrAdaptCoeff((traceFreqData).', (trackedVal).', freqExtr);
            
        else
            % -- Use matlab
            method = 'cubic'; % ; 'nearest', 'linear', 'spline', 'pchip', 'cubic'
            [ValExtr] = interp1(traceFreqData,trackedVal,freqExtr,method,'extrap');
        end
        
        plot(freqExtr,ValExtr,'sk:','LineWidth',2.5,'MarkerSize',12);%,'MarkerFaceColor','k');

    end%if
        
    