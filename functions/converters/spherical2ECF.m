function r_ECEF = spherical2ECF(lat, lon, alt, R_E)
    % Radius from Earth's center
    r = R_E + alt;
    
    % Position in ECEF frame at t = 0
    x = r * cos(lat) * cos(lon);
    y = r * cos(lat) * sin(lon);
    z = r * sin(lat);

    r_ECEF = [x; y; z];

end