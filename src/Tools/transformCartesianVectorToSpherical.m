function vectorAspherical = transformCartesianVectorToSpherical(vectorAcartesian, x, y, z)
    %vectorAspherical
    %   Usage:
    %           vectorAspherical = transformCartesianVectorToSpherical(vectorAcartesian, x, y, z)
    %
    %   Input Arguments:
    %       vectorAcartesian
    %           The x, y and z components of the spatial vector A = [Ax, Ay, Az]
    %       x,y,z
    %           The cartesian point where the vector is defined
    %
    %   Output Arguments:
    %       vectorAspherical
    %           The spherical components of the vector A = [Ar Atheta Aphi]
    %
    %   Description:
    %       Transforms A vector given in Cartesian Coordinates, to one
    %       in Spherical Coordinates
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.04
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za
    %   =======================

    r = sqrt(x^2 + y^2 + z^2);
    theta = acos(x/r);
    phi = atan(y/x); 
    if(y == 0)
        phi = 0;
    elseif(x == 0)
        phi = pi/2*y/abs(y);    
    end

    Ax = vectorAcartesian(1);
    Ay = vectorAcartesian(2);
    Az = vectorAcartesian(3);

    Ar = Ax*sin(theta)*cos(phi) + Ay*sin(theta)*sin(phi) + Az*cos(theta);
    Atheta = Ax*cos(theta)*cos(phi) + Ay*cos(theta)*sin(phi) - Az*sin(theta);
    Aphi = -Ax*sin(phi) + Ay*cos(phi);

    vectorAspherical = [Ar Atheta Aphi];