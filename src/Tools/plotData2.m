function plotData2(Const, x1,y1,x2,y2, leg_tr1,leg_tr2,xlab,ylab,titleString,imgString)
    %plotData v1.0
    %   Date: 13.06.2013
    %   Usage:
    %       plotData2(x1,y1,x2,y2, leg_tr1,leg_tr2,xlab,ylab,titleString,imgString)
    %
    %   Input Arguments:
    %       Const
    %           Struct containing information about the directory structure
    %       x1,y1
    %           The X and Y data of trace 1 
    %       x2,y2
    %           The X and Y data of trace 2
    %       leg_tr1
    %           The legend entry for trace 1
    %       leg_tr2
    %           The legend entry for trace 2
    %       xlab
    %           The label for the x-axis
    %       ylab
    %           The label for the y-axis
    %       titleString
    %           The figure title   
    %       imgString
    %           The name of the *.png and *.fig file of the image - will be
    %           stored in Const.OutputDirName
    %   Output:
    %       *.fig and *.png image files
    %
    %   Description:
    %       Plots two data sets / traces (x1,y1) and (x2,y2) on an
    %       graph and saves a *.fig and *.png image
    %
    %   =======================
    %   Written by Danie Ludick on June 13, 2013
    %   Last updated on June 21, 2013.
    %   EM Systems & Software (Pty) Ltd.
    %   Email: dludick.emss.co.za
    
    error(nargchk(11,11,nargin));

    fig = figure;
    hold on;
    grid on;
    box on;
    plot(x1,y1,'-.xb','LineWidth',4.0);
    plot(x2,y2,'-ok','LineWidt',4.0);
    leg = legend(leg_tr1,leg_tr2,2);
    %set(leg,'location','northwest'); %Not working
    xlabel(xlab);
    ylabel(ylab);
    title(titleString);
    
    % Save a *.png in the output directory
    pngName = strcat('/',imgString,'.png');
    saveas(fig,strcat(Const.OutputDirName,pngName),'png');
    
    % Save a *.fig in the output directory
    figName = strcat('/',imgString,'.fig');
    saveas(fig,strcat(Const.OutputDirName,figName),'fig');
    