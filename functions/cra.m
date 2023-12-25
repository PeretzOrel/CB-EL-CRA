function r = cra(params)
    x = []; y =[]; z = [];
    for i = 1:length(params.Rm)
        theta_ring = 2*pi*(0:params.Nm(i)-1)/params.Nm(i);
        x = [x; params.Rm(i)*cos(theta_ring)'];
        y = [y; params.Rm(i)*sin(theta_ring)'];
        z = [z; zeros(params.Nm(i),1)];
    end
    r = [x, y, z];
end