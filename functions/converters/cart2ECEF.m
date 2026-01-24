function [xp,yp,zp,vxp,vyp,vzp] = cart2ECEF(x,y,z,vx,vy,vz,t,mu,omega_E)
    % Latitude and longitude calculations
    a = zeros(1,length(t));
    ecc = a;
    inc = a;
    raan = a;
    argp = a;
    f = a;
    for j=1:length(t)
        [a(j),ecc(j),inc(j),raan(j),argp(j),f(j)] = cart2kep(x(j),y(j),z(j),vx(j),vy(j),vz(j),mu);
    end

    lag = omega_E*t;
    
    xp = a;
    yp = a;
    zp = a;
    vxp = a;
    vyp = a;
    vzp = a;
    for j=1:length(t)
        [xp(j),yp(j),zp(j),vxp(j),vyp(j),vzp(j)] = kep2cart(a(j),ecc(j),inc(j),raan(j)-lag(j),argp(j),f(j),mu);
    end
end