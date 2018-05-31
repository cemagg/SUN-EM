function plotData1(Const, x1,y1, xlab,ylab, titleString,imgString, fcdString)
    %plotData v2.0
    %   Date: 2016-06-22
    %   Usage:
    %       plotData1(x1,y1, xlab,ylab, titleString,imgString)
    %
    %   Input Arguments:
    %       Const
    %           Struct containing information about the directory structure
    %       x1,y1
    %           The X and Y data of trace 1
    %       xlab
    %           The label for the x-axis
    %       ylab
    %           The label for the y-axis
    %       titleString
    %           The figure title
    %       imgString
    %           The name of the *.png and *.fig file of the image - will be
    %           stored in Const.OutputDirName
    %       fcdString
    %           The name of the *.fcd filename.
    %   Output:
    %       *.fig and *.png image files. Also now add the option to store the plot
    %       data in a FEKO *.fcd file format (feko-connect data) file format that can
    %       be imported and plotted in POSTFEKO.
    %
    %   Description:
    %       Plots one data set / trace (x1,y1) on a
    %       graph and saves a *.fig and *.png image
    %
    %   =======================
    %   Written by Danie Ludick on June 13, 2013
    %   Last updated on June 21, 2013.
    %   EM Systems & Software (Pty) Ltd.
    %   Email: dludick.emss.co.za

    error(nargchk(8,8,nargin));

    figure;
    hold on;
    grid on;
    box on;
    xlabel(xlab);
    ylabel(ylab);
    title(titleString);

    %fig = figure;    
    if (Const.plotSemiLogY)
        fig = semilogy(x1,y1,'-.xb','LineWidth',2.0);
    else
        fig = plot(x1,y1,'-.xb','LineWidth',2.0);
    end%if

    % The following is only required when using Octave.
    if (Const.is_octave)
        graphics_toolkit gnuplot 
    end

    % Save a *.png in the output directory
    pngName = strcat('/',imgString,'.png');
    saveas(fig,strcat(Const.OutputDirName,pngName),'png');

    % Save a *.fig in the output directory
    figName = strcat('/',imgString,'.fig');
    saveas(fig,strcat(Const.OutputDirName,figName),'fig');

    % Option to store this trace data to a FEKO-Connect Data file (for later processing in POSTFEKO)
    storeFCDdata(Const, x1,y1, fcdString)
