function r_ECEF = ECI2ECEF(r_ECI, t, omega_E)
% ECI2ECEF - Convert inertial position to Earth-fixed at time t
% Applies a rotation about Z by +theta in frame sense, which for vector
% transformation corresponds to Rz(-theta) acting on the vector.
% Ensures output has same shape (row/column) as input.

    theta = omega_E * t;
    Rz_minus = [cos(theta), -sin(theta), 0;
                sin(theta),  cos(theta), 0;
                0         ,  0         , 1];

    if isrow(r_ECI)
        r_ECEF = (Rz_minus * r_ECI.').';
    else
        r_ECEF = Rz_minus * r_ECI;
    end
end