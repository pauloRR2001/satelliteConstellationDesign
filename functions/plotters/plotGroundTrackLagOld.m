function plotGroundTrackLag(x,y,z,vx,vy,vz,t,vernalEq,mu,omega_E)
    % r_vec_xyz: nx3 matrix [x, y, z] coordinates (km)
   
    [x,y,z,~,~,~] = cart2ECEF(x,y,z,vx,vy,vz,t,mu,omega_E);

    r = sqrt(x.^2 + y.^2 + z.^2);

    phi = acos(z ./ r); % spherecial coord 
    lat = -rad2deg(phi - pi / 2); % Latitude, deg
    lon = rad2deg(atan2(real(y), real(x))); % Longitude, deg

    % compensate location of vernal equinox
    lon = lon+rad2deg(vernalEq);
    
    opts = struct();
    opts.Color = [0, 0.4470, 0.7410]; % Default line color
    opts.LineStyle = '-'; % Default line style
    opts.LineWidth = 1.5; % Default line width

    % Plot the background (totally stole this)
    load('topo.mat', 'topo'); % Built-in MATLAB topographic data
    topoplot = [topo(:,181:360), topo(:,1:180)]; % longitude range
    %figure;
    contour(-180:179, -90:89, topoplot, [0, 0], 'black'); % Plot coastlines
    hold on;

    % Handle longitude discontinuity
    discontinuities = find(abs(diff(lon)) > 180);
    lon(discontinuities + 1) = NaN; % Insert NaN to avoid jumps

    % Plot ground track
    plot(lon, lat, 'Color', opts.Color, 'LineStyle', opts.LineStyle, 'LineWidth', opts.LineWidth);
    plot(lon(1), lat(1), 'bo');
    plot(lon(end), lat(end), 'rx');

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
