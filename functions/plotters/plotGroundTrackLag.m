function plotGroundTrackLag(a,ecc,inc,raan,argp,t,vernalEq,mu,omega_E)
    
    n = sqrt(mu/a^3);
    M_vec = n*t;
    
    for j=1:length(t)
        [x_ECI,y_ECI,z_ECI,~,~,~] =...
            kep2cart(a,ecc,inc,raan,argp,M_vec(j),mu);

        r_ECEF(:,j) = ECI2ECEF([x_ECI,y_ECI,z_ECI],t(j),omega_E);
    end

    x_ECEF = r_ECEF(1,:);
    y_ECEF = r_ECEF(2,:);
    z_ECEF = r_ECEF(3,:);

    r = sqrt(x_ECEF.^2 + y_ECEF.^2 + z_ECEF.^2);

    phi = acos(z_ECEF ./ r); % spherecial coord 
    lat = -rad2deg(phi - pi / 2); % Latitude, deg
    lon = rad2deg(atan2(y_ECEF, x_ECEF)); % Longitude, deg

    % compensate location of vernal equinox
    lon = lon+rad2deg(vernalEq);
    
    %opts = struct();
    %opts.Color = [0, 0.4470, 0.7410]; % Default line color
    %opts.LineStyle = '-'; % Default line style
    %opts.LineWidth = 1.5; % Default line width

    % Plot the background (totally stole this)
    load('topo.mat', 'topo'); % Built-in MATLAB topographic data
    topoplot = [topo(:,181:360), topo(:,1:180)]; % longitude range
    %figure;
    contour(-180:179, -90:89, topoplot, [0, 0], 'black'); % Plot coastlines
    hold on;

    % % Handle longitude discontinuity
    discontinuities = find(abs(diff(lon)) > 180);
    lon(discontinuities + 1) = NaN; % Insert NaN to avoid jumps

    % Plot ground track
    plot(lon, lat, 'LineWidth', 1.2);
    %plot(lon(1), lat(1), 'bo');
    %plot(lon(end), lat(end), 'rx');

    axis equal;
    grid on;
    xlim([-180, 180]);
    xticks(-180:30:180);
    ylim([-90, 90]);
    yticks(-90:30:90);
    xlabel('Longitude [deg]', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Latitude [deg]', 'Interpreter', 'latex', 'FontSize', 14);
    title('Satellite Ground Track', 'Interpreter', 'latex', 'FontSize', 16);
end
