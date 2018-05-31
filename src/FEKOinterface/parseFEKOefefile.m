function [FEKO_Efield] = parseFEKOefefile(Const, efe_filename)
    %parseFEKOefefile
    %   Date: 2018.05.24
    %   Usage:
    %       [Const, FEKO_data] = parseFEKOefefile(Const)
    %
    %   Input Arguments:
    %       Const: A global struct containing general program flow settings.
    %       efe_filename:
    %              The particular *.efe filename that will be used    
    %
    %   Output Arguments:
    %       FEKO_Efield:
    %           Struct containing the E-field sampling points and values.
    %
    %   Description:
    %       Reads in a FEKO *.efe file and extracts the electric field values,
    %       as well as sampling points.
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.24
    %   Stellenbosch University
    %   Email: dludick.sun.ac.za

    narginchk(2,2);

    fid = fopen(efe_filename,'r');

    if fid == -1
        message_fc(Const,sprintf('Error reading FEKO *.efe file: %s',efe_filename));
        error(['Error reading FEKO *.efe file: %s' efe_filename]);
    end

    message_fc(Const,' ');
    message_fc(Const,...
        '------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Parsing the FEKO *.efe file'));
    message_fc(Const,sprintf('  *.efe file: %s',efe_filename));

    % Initialise the return values.
    FEKO_Efield = [];

    % ========================
    % Read the file type
    % ========================
    line=fgetl(fid);
    efe_line_data = strsplit(line);
    % Make sure we are working with the correct data (Electric field)
    if (~strcmp(efe_line_data{3},'Electric'))
        message_fc(Const,sprintf('Expecting electric field values'));
        error(['Expecting electric field values']);
    end%if

    % ========================
    % Read the file version
    % ========================
    line=fgetl(fid);
    efe_line_data = strsplit(line);
    % Make sure we are working with the correct data (Electric field)
    if (str2num(efe_line_data{3}) ~= Const.FEKO_efe_file_format)
        message_fc(Const,sprintf('Unsupported file format for *.efe file'));
        error(['Unsupported file format for *.efe file']);
    end%if
    
    % Skip a few lines until we see the "#Frequency line" below.
    frequency_line_found = false;
    while(~frequency_line_found)
        line=fgetl(fid);
        efe_line_data = strsplit(line);
        if(strcmp(efe_line_data{1},'#Frequency:'))
            frequency_line_found = true;
        end%if
    end%while (~frequency_line_found)
    
    % ========================
    % Read the Frequency
    % ========================
    % We have already split the frequency line into tokens
    %line=fgetl(fid);
    %efe_line_data = strsplit(line);
    FEKO_Efield.frequency = str2double(efe_line_data{2});

    % ========================
    % Read the co-ordinate system
    % ========================
    line=fgetl(fid);
    efe_line_data = strsplit(line);
    FEKO_Efield.coordinate_system = efe_line_data{3};

    % ========================
    % Read the number of R samples
    % ========================
    line=fgetl(fid);
    efe_line_data = strsplit(line);
    FEKO_Efield.number_r_samples = str2num(efe_line_data{5});

    % ========================
    % Read the number of Theta samples
    % ========================
    line=fgetl(fid);
    efe_line_data = strsplit(line);
    FEKO_Efield.number_theta_samples = str2num(efe_line_data{5});

    % ========================
    % Read the number of Phi samples
    % ========================
    line=fgetl(fid);
    efe_line_data = strsplit(line);
    FEKO_Efield.number_phi_samples = str2num(efe_line_data{5});

    % Skip a few lines again
    line=fgetl(fid);
    line=fgetl(fid);
    line=fgetl(fid);
    
    number_field_points = FEKO_Efield.number_r_samples*FEKO_Efield.number_theta_samples*FEKO_Efield.number_phi_samples;

    % Initialise the field-point values
    % Note: The frequency axis has not yet been defined.

    % Spherical co-ordinate (r, theta, phi)
    FEKO_Efield.r_samples_m = zeros(number_field_points,1);
    FEKO_Efield.theta_samples_deg = zeros(number_field_points,1);
    FEKO_Efield.phi_samples_deg = zeros(number_field_points,1);
        
    % E-field value (Er, Etheta, Ephi)
    FEKO_Efield.Er = complex(zeros(number_field_points,1));
    FEKO_Efield.Etheta = complex(zeros(number_field_points,1));
    FEKO_Efield.Ephi = complex(zeros(number_field_points,1));

    for efield_indx = 1:number_field_points
        line=fgetl(fid);
        efe_line_data = strsplit(line);
        
        FEKO_Efield.r_samples_m(efield_indx,1) = str2double(efe_line_data{2});
        FEKO_Efield.theta_samples_deg(efield_indx,1) = str2double(efe_line_data{3});
        FEKO_Efield.phi_samples_deg(efield_indx,1) = str2double(efe_line_data{4});
        
        FEKO_Efield.Er(efield_indx,1) = str2double(efe_line_data{5}) + 1i*str2double(efe_line_data{6});
        FEKO_Efield.Etheta(efield_indx,1) = str2double(efe_line_data{7}) + 1i*str2double(efe_line_data{8});
        FEKO_Efield.Ephi(efield_indx,1) = str2double(efe_line_data{9}) + 1i*str2double(efe_line_data{10});

    end%while end_flag == 0