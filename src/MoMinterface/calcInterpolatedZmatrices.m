function [Interpolated_Data] = calcInterpolatedZmatrices(Const, Solver_setup, zMatrices)
    %runEMsolvers
    %   Usage:
    %       [Interpolated_Data] = calcInterpolatedZmatrices(Const, Solver_setup, zMatricesFEKO)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containging settings of which solver to run
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details
    %       zMatrices
    %           The Z-matrices data (all the values from e.g. FEKO)
    %   Output Arguments:
    %      Interpolated_Data
    %           The Z-matrices data (that includes both the calculated
    %           values from e.g. FEKO, as well as interpolated values
    %
    %   Description:
    %       Reads in impedance Z-matrices over a frequency range. For a
    %       givem Z(m,n) element it then extracts a number of calculated
    %       points and determines the rest using interpolation.

    narginchk(3,3);

    % Initialise the return value to be the same as that of FEKO (to get
    % the general structure correctly setup.
    Interpolated_Data = []; 
    zMatricesINTERP  = zMatrices.values;


    %zMatricesINTERP = [];

    frequency = Solver_setup.frequencies.samples;
    numFreq = Solver_setup.frequencies.freq_num;

     % for freq = 1:numFreq 
    RWGmBasis = 200; %Const.numMoMbasis;
    RWGnBasis = 200; %Const.numMoMbasis;
    fstep = 2;
 

%empty the matrices and leave the diagonals, retain selected frequencies,
%50% retained
for freq = 1:fstep:numFreq
    for m = 1:RWGmBasis   
        for n = 1:RWGnBasis
            if m ~= n
            zMatricesINTERP(n,m,freq) = 0;
           end
        end 
    end
end





% % Improve the matrices[Zmn] for selected frequencies
for freq = 2:fstep:numFreq
    for m = 1:RWGmBasis   
        for n = 1:RWGnBasis
            if m ~= n
                
                FrequencySample = frequency(freq);
                lambda = physconst('LightSpeed')./FrequencySample;

                edge_m_X = Solver_setup.rwg_basis_functions_shared_edge_centre(m,1);
                edge_m_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(m,2);

                edge_n_X = Solver_setup.rwg_basis_functions_shared_edge_centre(n,1);
                edge_n_Y = Solver_setup.rwg_basis_functions_shared_edge_centre(n,2);

                Rmn = sqrt((edge_m_X - edge_n_X)^2 + (edge_m_Y - edge_n_Y)^2);

                zMatricesINTERP(n,m,freq) = zMatricesINTERP(n,m,freq)./exp(-1i*((2*pi)./lambda')*Rmn);
           end
        end 
    end
end

Interpolated_Data.Interpolate_Zmn = InterpolateZmn(Const, Solver_setup, zMatricesINTERP, zMatrices);



%  zMatricesFEKOcalulatedd = [];

%     
% % Zmn = zMatrices.values(1,200,1:5);
% for  m = 1:RWGmBasis
%   for  n= 1:RWGnBasis
%         
%     if m ~= n
% 
%           frequency = Solver_setup.frequencies.samples;
%           numFreq = Solver_setup.frequencies.freq_num;
%           FewFreqsamples = frequency(1:fstep:numFreq);
%           lambda = physconst('LightSpeed')./FewFreqsamples;
%           
%           zMatricesFEKO = zMatrices.values(m,n,1:fstep:numFreq);   % build 3D array of all of individuals to manipulate as one
%           zMatricesFEKO = reshape(permute(zMatricesFEKO,[3,2,1]),50,[]); % reshape vector to matrix
% 
%           zMatricesINTERP = [zMatricesINTERP; zMatricesFEKO];
% 
%          %Calculated zMatrices at selected frequencies
%          %50 evenly spaced frequencies
% 
%         % zMatricesFEKOcalulated = zMatrices.values(m,n,1:fstep:numFreq);
%          %zMatricesFEKOcalulatedd = [zMatricesFEKOcalulatedd; zMatricesFEKOcalulated];
%          
%          
%          %call the function for interpolation
%     end
%    end 
% end
% 
% %zMatricesInterpolated = (InterpolateZmn(Const, Solver_setup, zMatricesINTERP));
% 
% %Zmom = (calcZmn(Const, zMatrices, freq, 1, 1, ObservRWGs, SourceRWGs));
% 
