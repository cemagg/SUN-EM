function storeFCDdata(Const, x1,y1, fcdString)
    %plotData v1.0
    %   Date: 2016-06-27
    %   Usage:
    %       plotData1(Const, x1,y1, fcdString)
    %
    %   Input Arguments:
    %       Const
    %           Struct containing information about the directory structure
    %       x1,y1
    %           The X and Y data of trace 1
    %       fcdString
    %           The name of the *.fcd filename.
    %   Output:
    %       Saves data to an *.fcd data file (FEKO Connect data file) - CSV delimeted file.
    %
    %   Description:
    %       See above.
    %
    %   =======================
    %   Written by Danie Ludick on June 27, 2016
    %   Email: dludick@sun.ac.za

    error(nargchk(4,4,nargin));

    % Option to store this trace data to a FEKO-Connect Data file (for later processing in POSTFEKO)
    if (Const.store_to_fcd_file)
        % Open the *.fcd in the project directory.
        fcdFileName = strcat('/',fcdString,'.fcd');
        fidOut = fopen(strcat(Const.OutputDirName,fcdFileName), 'w');
        % Loop over the data and write it in CSV file-format
        for n=1:max(size(x1,1),size(x1,2))
            fprintf (fidOut, '%20.10e,%20.10e\n', x1(n), y1(n));
        end%for
        % Close the file
        fclose(fidOut);
    end%if
