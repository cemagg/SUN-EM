function [segments] = parseFEKOoutfileSegments(Const)
    %parseFEKOoutfileSegments v1.0
    %   Date: 12.06.2013
    %   Usage:
    %       [segments] = parseFEKOoutfileSegments(Const)
    %
    %   Input Arguments:
    %       Const: A global struct containing:
    %       FEKOoutfilename
    %           FEKO *.out filename (e.g. 'yagi.out')
    %
    %   Output Arguments:
    %       Segments:
    %           The struct containing the segment data
    %
    %   Description:
    %       Reads in a FEKO *.out file which (potentially) contains
    %       segments and extracts the segment data (start and end
    %       positions & radius) and stores it for further processing.
    %
    %   =======================
    %   Written by Danie Ludick on June 12, 2013
    %   Last updated on June 12, 2013.
    %   EM Systems & Software (Pty) Ltd.
    %   Email: dludick.emss.co.za

    %   Please note that additional information on reading the *.out files using FEKO can be found 
    %   on the FEKO website:
    %   http://www.feko.info/
    
    error(nargchk(1,1,nargin));

    fid = fopen(Const.FEKOoutfilename,'r');

    if fid == -1
        message(Const,sprintf('Error reading FEKO *.out file: %s',Const.FEKOoutfilename));
        error(['Error reading FEKO *.out file: %s' Const.FEKOoutfilename]);
    end

    message(Const,sprintf('Reading segment data from *.out file: %s',Const.FEKOoutfilename) )
    
    end_flag = 0;
    end_seg_read = 0;
    seg_num = 0;
    impedance_real = [];
    
    segments.start_xcordinates = [];
    segments.start_ycordinates = [];
    segments.start_zcordinates = [];
    segments.end_xcordinates = [];
    segments.end_ycordinates = [];
    segments.end_zcordinates = [];
    segments.length = [];
    segments.radius = [];
    
    while end_flag == 0
        line=fgetl(fid);
        
        % Check for end of file:
        if strcmp(line,'                    SUMMARY OF REQUIRED TIMES IN SECONDS');
            end_flag = 1;
        end%if
        
        % Start processing the segment data:
        if strcmp(line,'                              DATA OF THE SEGMENTS');
            line=fgetl(fid);
            line=fgetl(fid);
            line=fgetl(fid);
            g = strfind(line,'medium');
            if g > 0
                % The next line is the start of the co-ordinates (two lines)
                while end_seg_read == 0
                    line=fgetl(fid);
                    if (length(line)==0)
                        end_seg_read =1;
                    else
                        seg_num = seg_num + 1;
                        % Read the start position (X,Y,Z in [m])         
                        segments.start_xcordinates(seg_num) = str2num(line(19:29));
                        segments.start_ycordinates(seg_num) = str2num(line(31:41));
                        segments.start_zcordinates(seg_num) = str2num(line(43:53));
                        % Read the end position, length and radius (all in [m])
                        line=fgetl(fid);
                        segments.end_xcordinates(seg_num) = str2num(line(19:29));
                        segments.end_ycordinates(seg_num) = str2num(line(31:41));
                        segments.end_zcordinates(seg_num) = str2num(line(43:53));
                        segments.length(seg_num) = str2num(line(57:66));
                        segments.radius(seg_num) = str2num(line(70:79));
                   end
                end%while end_seg_read == 0 
            end%if
        end
        
    end%while end_flag == 0
    fclose('all');
    
    segments.seg_num = seg_num;
    
    message(Const,sprintf('Read: %d segments',seg_num));
    message(Const,sprintf('Finished processing the *.out file: %s',Const.FEKOoutfilename));
