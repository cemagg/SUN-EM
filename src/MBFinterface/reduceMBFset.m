function JReduced = reduceMBFset(Const, JFull, fcdString)
    %reduceMBFset
    %   Date: 2016-06-27
    %   Usage:
    %       reduceMBFset(Const, JFull, fcdString)
    %
    %   Input Arguments:
    %       Const
    %           Struct containing information about the directory structure
    %       JFull
    %           The column augmented data on which the SVD will be applied.
    %       fcdString
    %           The name of the *.fcd filename.
    %   Output:
    %       Orthonormalised a set of MBFs by applying the SVD and retaining only
    %       the number after applying a preset threshold.
    %
    %   Description:
    %       See above.
    %
    %   =======================
    %   Written by Danie Ludick on June 27, 2016
    %   Email: dludick@sun.ac.za
    %
    %   Contributions:
    %   22-Oct-2010 Rob Maaskant, ASTRON - original
    %   18-Aug-2015 Danie Ludick, Stellenbosch University (modified)

    LOCAL_DEBUG = false;

    % TO-DO: Danie, obtain a better understanding for the computational costs of the SVD in Matlab
    [A,B,C] = svd(JFull,0);

    Tol = Const.MBFthreshold;

    % 2016-06-22: Added now a flag here to retain all the orthonormalised bases, not just as specified 
    % by the threshold.
    if (Tol == -1)
        keep_all_bases = true;
    else
        keep_all_bases = false;    
    end%if

    if size(B,1)>1

        if (keep_all_bases)
            NumberOfRetainedBases=size(B,1);
        else
            NumberOfRetainedBases=sum(diag(B)./max(diag(B))>(1/Tol),1);
        end%if
        JReduced=A(:,1:NumberOfRetainedBases);

        if LOCAL_DEBUG
            figure
            diagB=diag(B);
            semilogy([1:length(diagB)],diagB./max(diagB),'LineWidth',2,'Color','b','Marker','s');
            % Also store the data to an FCD data-file:
            storeFCDdata(Const, [1:length(diagB)],diagB./max(diagB), fcdString);
            storeFCDdata(Const, [1:length(diagB)],(1/Tol) .* ones(size(diagB)), strcat('threshold_',fcdString));
            hold on
            if (~keep_all_bases)
                % Add also a threshold here if we are not going to keep all the basis functions.
                line([1 length(diagB)],[1/Tol 1/Tol],'LineStyle','--','LineWidth',2,'Color','r');
            end%if
            %xlim([1 100]);
            %xlim([1 size(JFull,2)])
            grid on
            xlabel('n    [-]','FontSize',12,'FontWeight','bold')
            ylabel('\sigma_{n}/\sigma_{max}    [-]','FontSize',12,'FontWeight','bold')
            if (~keep_all_bases)
                % Only add the legend entries, if we have specified a threshold to limit the MBFs.
                legend('Singular Values','Threshold')
            end%if
            set(gca,'FontSize',12)
            set(gca,'FontWeight','bold')
        end
    else
        NumberOfRetainedBases=1;
        JReduced=A(:,1:NumberOfRetainedBases);
    end