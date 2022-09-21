   function [zMatrices] = readFEKOZMatrixFromFile(Const, mat_file_name)
    %readFEKOZMatrixFromFile
    %   Usage:
    %       [zMatrices] = readFEKOZMatrixFromFile(Const)
    %
    %   Input Arguments:
    %       Const
    %           A global struct with some program flow parameters and constants
    %       mat_file_name
    %           FEKO matrix filename (e.g. 'yagi.mat')
    %
    %   Output Arguments:
    %       zMatrices, a struct containing:
    %       values
    %           Set of P FEKO Z matrices (MxNxP), representing an MxN
    %           matrix at P frequencies. (The file format is setup such
    %           that the Z matrix may be rectangular instead of square, but
    %           this is highly unlikely.)
    %       fileVersion
    %           The file version of the *.str file
    %       md5Check
    %           The MD5 checksum of the *.str file
    %       filePrecision
    %           The precision file, i.e. 'single' or 'double'
    %       mBasis
    %           The number of rows of the MxN (MoM) matrix
    %       nBasis
    %           The number of collumns of the MxN (MoM) matrix
    %       numFreq
    %           The number of frequencies, i.e. Z-mat blocks read
    %
    %   Description:
    %       Reads in a FEKO binary Z matrix file which contains either one
    %       or more matrices, each at a different frequency. Note, FEM is
    %       not supported here.
    %
    %       This code is based on sample Fortran code given by FEKO Support
    %       USA (specifically, Dr. Rensheng (Ray) Sun, ray@emssusa.com) on
    %       June 10, 2008.
    %
    %   =======================
    %   Originally written by Bryan Raines on June 11, 2008
    %   Last updated on June 11, 2008.
    %   ElectroScience Laboratory at The Ohio State University
    %   Email: rainesb@ece.osu.edu

    %   =======================
    %   Adapted by Danie Ludick on 26th of July 2012
    %   to account also for the fact that the *.mat file can be generated
    %   with 32-bit or 64-bit FEKO
    %   Last updated on July 26, 2012.
    %   EMSS-SA (Pty) Ltd
    %   Email: dludick@emss.co.za
    %   =======================
    %
    %   Adapted by Evan Lezar on 5th of December 2013
    %   to automatically handle the various padding options of FEKO mat files.
    %   Also corrected a bug in the reading of single precision data.
    %   Last updated on December 5, 2013.
    %   EMSS-SA (Pty) Ltd
    %   Email: lezar@emss.co.za

    %   Set constants related to which FEKO executables was used to generate the
    %   *.mat file, i.e. either 32-bit and 64-bit FEKO executables. Please
    %   note that additional information on reading the *.mat files using
    %   FEKO can be found on the FEKO website:
    %   http://www.feko.info/support/helpcenter/how-to/how-to-read-the-.mat-.lud-.rhs-files-and-.str-files

    %   (Note, as of FEKO Suite 6.2 the BYTE_PADDING for the *.mat file using both 32-bit and
    %   64-bit FEKO exectubles will be 4 Bytes) and the following settings
    %   should be fine for both 32-bit and 64-bit.
    PADDING_TYPE = 'auto';
    % If the automatic determination does not work as expected, the following parameters can be tried.
    % PADDING_TYPE = 'int32' % For 32-bit FEKO executables
    % PADDING_TYPE = 'int64' % For 64-bit FEKO executables (6.1.1 and earlier)

    narginchk(2,2);

    % After the *.mat file name is set above, open it now.
    fid = fopen(mat_file_name,'r','l');

    if fid == -1
        error(['Error reading FEKO matrix file: ' mat_file_name]);
        mesage(Const,sprintf('Error read: %s',mat_file_name));
    end

    if PADDING_TYPE == 'auto'
        auto_used = true;
        % First try the padding type of int64 to support 6.1.1 and earlier.
        PADDING_TYPE = 'int64';
    else
        auto_used = false;
    end

    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Reading *.mat file from: %s',mat_file_name));

    % Attempt to read the file version
    %File version
    while true
        reclen = fread(fid,1,PADDING_TYPE);
        zMatrices.fileVersion = fread(fid,1,'int32');
        reclen_end = fread(fid,1,PADDING_TYPE);
        % The start and end record length should match.
        if (reclen ~= reclen_end)
            if (auto_used && strcmp(PADDING_TYPE,'int64'))
                % If we have just tried padding type 'int64', then we must also try
                % the padding type of 'int32'
                PADDING_TYPE = 'int32';
                % Reset the file to read the value again.
                fseek(fid,0,'bof');
            else
                % This is an error state.
                message_fc(Const,'Invalid Record Length');
                error ('Invalid Record Length');
            end
        else
            break
        end
    end

    %MD5 Check
    reclen = fread(fid,1,PADDING_TYPE);
    md5Check = char(fread(fid,reclen,'uchar'));
    zMatrices.md5Check = '';
    for i=1:size(md5Check,1)
        zMatrices.md5Check = strcat(zMatrices.md5Check, md5Check(i));
    end
    reclen_end = fread(fid,1,PADDING_TYPE);
    % The start and end record length should match.
    if (reclen ~= reclen_end)
        message_fc(Const,'Invalid Record Length');
        error ('Invalid Record Length');
    end

    %For version >= 2, read file precision field
    file_single_precision = false;
    if version >= 2
        reclen = fread(fid,1,PADDING_TYPE);
        filePrecision = fread(fid,1,'int32');
        reclen_end = fread(fid,1,PADDING_TYPE);
        % The start and end record length should match.
        if (reclen ~= reclen_end)
            message_fc(Const,'Invalid Record Length');
            error ('Invalid Record Length');
        end
        % Check to see if the file was written in single precision.
        % This is written out as a Fortran LOGICAL
        % therefore true has the value 0xFFFFFFFF
        %      and false has the value 0x00000000
        if (filePrecision ~= 0)
            file_single_precision = true;
        end
    end

    if (file_single_precision)
        zMatrices.filePrecision = 'single';
        precisionString = 'single';
        element_size = 4;
    else
        zMatrices.filePrecision = 'double';
        precisionString = 'double';
        element_size = 8;
    end%if

    %Read number of elements (MxN matrix)
    %(1) Read M
    reclen = fread(fid,1,PADDING_TYPE);
    MM = fread(fid,1,'int32');
    reclen_end = fread(fid,1,PADDING_TYPE);
    % The start and end record length should match.
    if (reclen ~= reclen_end)
        message_fc(Const,'Invalid Record Length');
        error ('Invalid Record Length');
    end
    zMatrices.mBasis = MM;

    %(2) Read N
    reclen = fread(fid,1,PADDING_TYPE);
    NN = fread(fid,1,'int32');
    reclen_end = fread(fid,1,PADDING_TYPE);
    % The start and end record length should match.
    if (reclen ~= reclen_end)
        message_fc(Const,'Invalid Record Length');
        error ('Invalid Record Length');
    end
    zMatrices.nBasis = NN;

    matrixSize = MM*NN;

    freqIndex = 1;
    zMatrices.values = [];
    rawData = zeros(2*MM,NN);
    while ~feof(fid)
        try
            for colIndex = 1:NN
                % read the data for the column
                reclen = fread(fid,1,PADDING_TYPE);
                rawData(:,colIndex) = fread(fid,reclen/element_size,precisionString);
                reclen_end = fread(fid,1,PADDING_TYPE);
                % The start and end record length should match.
                if (reclen ~= reclen_end)
                    message_fc(Const,'Invalid Record Length');
                    error ('Invalid Record Length');
                end
            end
        catch
            % In this case an error has probably occurred.
            if freqIndex == 1
                message_fc(Const,'No matrix data read');
                error (['No matrix data read']);
            end
            break;
        end

        zMatrices.values(1:MM,1:NN,freqIndex) = complex(rawData(1:2:2*MM,:),rawData(2:2:2*MM,:));
        freqIndex = freqIndex+1;
    end

    fclose(fid);

    % Set the number of frequency samples that we have calculated here.
    zMatrices.numFreq = freqIndex-1;

    message_fc(Const,sprintf('Read Zmn file version %d',zMatrices.fileVersion));
    message_fc(Const,sprintf('Dim. of Zmn is m=%d and n=%d',zMatrices.mBasis,zMatrices.nBasis));
    message_fc(Const,sprintf('Zmn is stored in %s precision',zMatrices.filePrecision));
    message_fc(Const,sprintf('Read: %d frequencies (number of Z-mats.)',zMatrices.numFreq));
    message_fc(Const,sprintf('Finished processing the *.mat file: %s',Const.FEKOmatfilename));