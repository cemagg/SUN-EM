function EfieldAtPointSpherical =  calculateEfieldAtPointRWG(Const, r, theta_degrees, phi_degrees, ...
    Solver_setup, Isol)
    %calculateEfieldAtPointRWG
    %   Usage:
    %       [EfieldAtPointSpherical] =  calculateEfieldAtPointRWG(r, theta, phi, nodes, ...
    %           sharedEdgesList, triangleDataList, Isol, k )
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       r, theta_degrees, phi_degrees
    %           Spherical co-ordinate of the observation point. Note the angles are in degrees.
    %       Solver_setup (struct)
    %           Struct containing the frequency range, nodal co-ordinates of the 
    %           triangles, as well as information on the shared edges, RWG basis functions, 
    %           etc. 
    %       Isol
    %           The solution vector, i.e. the MoM expansion coefficients
    %   Output Arguments:
    %       EfieldAtPointSpherical
    %           The E-field strength calculated at the point (r,theta,phi)
    %           EfieldAtPoint = A r^ + B theta^ + C phi^
    %
    %   Description:
    %       Calculates the E-field value at a certain point
    %       (r,theta,phi). The E-field value is based on the dipole
    %       model, as discussed in [1] and [2].
    %
    %   TO-DO: Z-displaced atennas is not working correctly
    %
    %   References:
    %   [1] S. Makarov, "MoM Antenna Simulations with Matlab: RWG Basis Functions"
    %   [2] Balanis, "Antenna Theory: Analysis and Design (3rd Edition)"
    %   =======================
    %   Written by Danie Ludick on 2018.05.04
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    % Initialisations
    % Note: for each frequency, we will have a different field
    EfieldAtPointCartesian = zeros(Solver_setup.frequencies.freq_num, 3);
    EfieldAtPointSpherical = zeros(Solver_setup.frequencies.freq_num, 3);

    % Transform the point P = (r,theta,phi) to P = (x,y,z)
    [Px,Py,Pz] = transformSpericalCoordinateToCartesian(r,theta_degrees*Const.DEG2RAD,phi_degrees*Const.DEG2RAD);

    for freq_index = 1:Solver_setup.frequencies.freq_num

        % Wavelength
        lambda = Const.C0/Solver_setup.frequencies.samples(freq_index);
        k = 2*pi/lambda;
        K = (k*Const.ETA_0/(4*pi));
        %K = 1.0;

        % The E-field at the point is calculated as the superposition of the individual 
        % RWG (dipole) contributions. See definitions for m, M and C (Eqs. 9a and 9b)
        for rwg_bf_index = 1:Solver_setup.num_metallic_edges

            % --------------------------------
            % calculate m (the dipole moment vector) associate with the mth edge: 
            %     m = lm*Im*[ (r_mc^-) - (r_mc^+)]
            % where "c" denotes the triangle midpoint.
            % --------------------------------
            lengthM = Solver_setup.rwg_basis_functions_length_m(rwg_bf_index);
            Im = Isol(rwg_bf_index);

            triangle_plus = Solver_setup.rwg_basis_functions_trianglePlus(rwg_bf_index);
            triangle_minus = Solver_setup.rwg_basis_functions_triangleMinus(rwg_bf_index);

            rMplusX  = Solver_setup.triangle_centre_point(triangle_plus,1);
            rMminusX = Solver_setup.triangle_centre_point(triangle_minus,1);
            
            rMplusY  = Solver_setup.triangle_centre_point(triangle_plus,2);
            rMminusY = Solver_setup.triangle_centre_point(triangle_minus,2);
            
            rMplusZ  = Solver_setup.triangle_centre_point(triangle_plus,3);
            rMminusZ = Solver_setup.triangle_centre_point(triangle_minus,3);

            rwgDipoleCentreXYZ  = lineCentre([rMplusX, rMplusY, rMplusZ],[rMminusX, rMminusY, rMminusZ]);
            rwgDipoleCentre = rwgDipoleCentreXYZ;
            
            rL = sqrt( (rwgDipoleCentre(1) - Px)^2 + (rwgDipoleCentre(2) - Py)^2 + ...
                (rwgDipoleCentre(3) - Pz)^2 );
            
            rLvecX = Px - rwgDipoleCentreXYZ(1);
            rLvecY = Py - rwgDipoleCentreXYZ(2);
            rLvecZ = Pz - rwgDipoleCentreXYZ(3);        
                                                            
            C = (1/(rL*rL)) * (1 + 1/(1j*k*rL));

            mX = lengthM*Im*(rMminusX - rMplusX);
            mY = lengthM*Im*(rMminusY - rMplusY);
            mZ = lengthM*Im*(rMminusZ - rMplusZ);   

            % --------------------------------
            % calculate M (Eq. 9a in [1])
            % --------------------------------
            mDotP = rLvecX*mX +  rLvecY*mY +  rLvecZ*mZ;

            % Mx = (1/(rL*rL)) * Px * (mDotP);
            % My = (1/(rL*rL)) * Py * (mDotP);
            % Mz = (1/(rL*rL)) * Pz * (mDotP);

            Mx = (1/(rL*rL))  * rLvecX * (mDotP);
            My = (1/(rL*rL))  * rLvecY * (mDotP);
            Mz = (1/(rL*rL))  * rLvecZ * (mDotP);

            % calculate the E-field (x,y,z) due to this RWG Basis Function. Note, the phase is referenced 
            % relative to the centre of the dipole. Account for that below.
            phaseTerm = exp(-1j*k*rL);  
            
            EfieldCurrentRWGelementX = K * ( (Mx - mX)*(1j*k*(1/rL) + C) + 2*Mx*C) * phaseTerm;
            EfieldCurrentRWGelementY = K * ( (My - mY)*(1j*k*(1/rL) + C) + 2*My*C) * phaseTerm;
            EfieldCurrentRWGelementZ = K * ( (Mz - mZ)*(1j*k*(1/rL) + C) + 2*Mz*C) * phaseTerm;

            EfieldAtPointCartesian(freq_index, 1) = EfieldAtPointCartesian(freq_index, 1) + EfieldCurrentRWGelementX;
            EfieldAtPointCartesian(freq_index, 2) = EfieldAtPointCartesian(freq_index, 2) + EfieldCurrentRWGelementY;
            EfieldAtPointCartesian(freq_index, 3) = EfieldAtPointCartesian(freq_index, 3) + EfieldCurrentRWGelementZ;
        end

        % Kdl - This is a constant 2*pi that has been added to get SUN-EM results the same as that of
        % FEKO. TO-DO: DL, check why this is the case. Peak vs. RMS? Cannot be, as that is a factor 1/2
        % or 1/sqrt(2).
        Kdl = 2*pi;

        EfieldAtPointCartesian = EfieldAtPointCartesian ./ Kdl;

        % transform the E-field vector, to one in Sperical Coordinates
        EfieldAtPointSpherical(freq_index,:) = transformCartesianVectorToSpherical(EfieldAtPointCartesian(freq_index,:), Px, Py, Pz);

    end%for freq_index = 1:Solver_setup.freq_num


