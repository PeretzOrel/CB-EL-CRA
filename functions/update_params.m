function params = update_params(params)
    params.lambda = params.c/max(params.fband);
    params.f_grid = linspace(params.fband(1), params.fband(2), params.K).';
    params.theta_BW_grid = linspace(0, params.theta_BW/2, 20).';
    params.theta_grid = linspace(-90, 90, 200).';
    params.M = length(params.Rm);
%     params.Nm = ceil(4*pi*params.Rm*max(params.fband)/params.c);
%     params.Nm(params.Nm == 0) = 1;
    params.N = sum(params.Nm); % num of elements
    params.r = cra(params);
    params.d = @(f, theta, phi) exp(-1j*(2*pi*f)*(1/params.c)*params.r*[sind(theta).*cosd(phi); sind(theta).*sind(phi); cosd(theta)]);
    params.d = @(f, theta) params.d(f, theta, 0);
    params.d_rings = @(f, theta) besselj(0, 2*pi.*params.Rm(:)*f/params.c.*sind(theta(:).'));
    params.T = select_rings_mat(1:params.M, params);
    params.T_normalized = params.T./sum(params.T);
end
