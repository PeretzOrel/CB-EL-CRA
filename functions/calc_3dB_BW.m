function bw = calc_3dB_BW(H, params)
    
    theta_pos = params.theta_grid(params.theta_grid >= 0);
    if params.K > 1
        bp = B(H, params.d, params.f_grid, theta_pos);
        bw = zeros(params.K, 1);
        for i = 1:params.K
            h = H(:, i);
            
            bpi = abs(bp(i, :)) - db2mag(-3);

            theta_h_idx = find(bpi <= 0, 1, 'first');
            theta_l_idx = theta_h_idx - 1;
            
            if isempty(theta_h_idx)
                bw(i) = 180;
                continue;
            end
            
            f = @(theta) abs(h'*params.d(params.f_grid(i), theta)) - db2mag(-3);
            p = bisection(f, theta_pos(theta_l_idx), theta_pos(theta_h_idx));
            bw(i) = 2*p;
        end
    else
        bp = B(H, params.d, params.f_grid, theta_pos);
        h = H;

        bpi = abs(bp) - db2mag(-3);

        theta_h_idx = find(bpi <= 0, 1, 'first');
        theta_l_idx = theta_h_idx - 1;

        if isempty(theta_h_idx)
            bw = 180;
        else
            f = @(theta) abs(h'*params.d(params.f_grid, theta)) - db2mag(-3);
            p = bisection(f, theta_pos(theta_l_idx), theta_pos(theta_h_idx));
            bw = 2*p;
        end
    end
end