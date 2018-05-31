function [feko_ffe_results] = parseFEKOffefile(Const, ffe_filename)
    %parse_FEKO_ffe_file
    %   Date: 2017.12.13
    %   Usage:
    %       [feko_ffe_results] = parse_FEKO_ffe_file(Const,ffe_filename)
    %
    %   Input Arguments:
    %       Const: A global struct containing:
    %       ffe_filename:
    %              The particular *.ffe filename that will be used
    %
    %   Output Arguments:
    %       feko_ffe_results:
    %           A struct containing the theta and phi samples, as well as the Etheta and Ephi components
    %           An N x 4 matrix containing for each of the N field-points the following:
    %               theta (in degrees), phi (in degrees), E_theta (complex), E_phi (complex)
    %
    %   Description:
    %       Reads in a FEKO *.ffe file (i.e. a far field) extracts the angles and field values
    %       in theta, phi, Etheta and Ephi
    %
    %   =======================
    %   Written by Danie Ludick on 2017.12.13
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    narginchk(2,2);

    % First check that the *.ffe file actually exists:
    if (~exist(ffe_filename,'file'))
        message_fc(Const,sprintf('Error locating *.ffe file: %s', ffe_filename));
        error(['Error locating *.ffe file: ' ffe_filename]);
    end%if

    fid = fopen(ffe_filename,'r');

    if fid == -1
        message_fc(Const,sprintf('Error reading FEKO *.ffe file: %s',ffe_filename));
        error(['Error reading FEKO *.ffe file: %s' ffe_filename]);
    end

    message_fc(Const,sprintf('Parsing the FEKO *.ffe file'));
    message_fc(Const,sprintf('  *.ffe file: %s',ffe_filename) );

    % Initialise the return struct
    feko_ffe_results = [];

    % ----------- Start reading the header information - extract the number of samples
    num_theta_phi_samples = 0;
    num_theta = 0;
    num_phi = 0;

    while (num_theta_phi_samples == 0)

        line=fgetl(fid);

        % TO-DO: Perhaps also other data later on ... e.g. 
        % file-version, frequency, etc. for now, only read theta and phi
        % samples

        % -------------------------------------------------
        % -- Parse the number of theta samples:
        g = strfind(line,'#No. of Theta Samples:'); 
        if (g > 0)
            theta_samples = strsplit(line);
            num_theta = str2num(theta_samples{5});
        end%if

        % -------------------------------------------------
        % -- Parse the number of phi samples:
        g = strfind(line,'#No. of Phi Samples:'); 
        if (g > 0)
            phi_samples = strsplit(line);
            num_phi = str2num(phi_samples{5});
        end%if

        num_theta_phi_samples = num_theta*num_phi;

    end%while end_flag == 0

    % move now to the start of the Etheta and Ephi samples
    %message_fc(Const,sprintf('              Reading now %d angles (theta = %d x phi = %d)', ...
    %    num_theta_phi_samples, num_theta, num_phi) );

    % Move past a couple of more lines in the header.
    line=fgetl(fid);
    line=fgetl(fid);
    line=fgetl(fid);

    % Store also the number of far field samples
    feko_ffe_results.number_theta_samples = num_theta;
    feko_ffe_results.number_phi_samples = num_phi;

    % Allocate some space for the return vectors (in the struct):
    feko_ffe_results.theta_samples_deg = zeros(num_theta_phi_samples,1);
    feko_ffe_results.phi_samples_deg = zeros(num_theta_phi_samples,1);
    feko_ffe_results.Etheta = complex(zeros(num_theta_phi_samples,1));
    feko_ffe_results.Ephi = complex(zeros(num_theta_phi_samples,1));
    
    % Read the theta, phi, real(Etheta), imag(Etheta), real(Ephi), imag(Ephi) data
    for sample = 1:num_theta_phi_samples
        line=fgetl(fid);
        data = strsplit(line);

        % Theta in degrees
        feko_ffe_results.theta_samples_deg(sample) = str2double(data{2});

        % Phi in degrees
        feko_ffe_results.phi_samples_deg(sample) = str2double(data{3});

        % Etheta
        feko_ffe_results.Etheta(sample) = str2double(data{4}) + 1i*str2double(data{5});

        % Ephi
        feko_ffe_results.Ephi(sample) = str2double(data{6}) + 1i*str2double(data{7});
    end%for

    % -- Finished, close the file again
    fclose(fid);

    return;