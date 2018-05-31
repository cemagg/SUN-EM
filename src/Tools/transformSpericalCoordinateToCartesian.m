function [x,y,z] = transformSpericalCoordinateToCartesian(r,theta_radians,phi_radians)
    %transformSpericalCoordinateToCartesian
    %   Usage:
    %           [x,y,z] = transformSpericalCoordinateToCartesian(r,theta,phi)
    %
    %   Input Arguments:
    %       r, theta_radians, phi_radians
    %           The spherical co-ordinate. Note, angle is in radians
    %
    %   Output Arguments:
    %       x,y,z
    %           The cartesian point
    %
    %   Description:
    %       Transforms Spherical Coordinates to Cartesian
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.04
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za
    %   =======================

    x = r*sin(theta_radians)*cos(phi_radians);
    y = r*sin(theta_radians)*sin(phi_radians);
    z = r*cos(theta_radians);
    