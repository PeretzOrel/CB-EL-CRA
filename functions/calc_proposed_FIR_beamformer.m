function [costval, optval] = calc_proposed_FIR_beamformer(params)
    G = [ones(params.K, 1), 2*cos(2*pi*params.f_grid*(1:params.L/2)*(1/params.Fs))].';
%     tic;
    cvx_begin quiet
    cvx_solver Mosek % SDPT3, sedumi , mosek
        variable B_opt(params.L/2+1, params.M);
        variable A(params.M, params.K);
        A_ext = params.T_normalized*A;

        dist = repmat(squareform(pdist(params.r)), 1, 1, params.K);
        ff(1, 1, :) = params.f_grid;
        ff = repmat(ff, params.N, params.N, 1);
        Gamma = sinc(2*dist.*ff/params.c);

        obj_fun = 0;
        for i = 1:params.K
            if params.alpha == 1
                obj_fun = obj_fun + quad_form(A_ext(:, i), Gamma(:, :, i));
            elseif params.alpha == 0
                obj_fun = obj_fun + norm(A_ext(:, i));
            else
                obj_fun = obj_fun + ...
                    params.alpha*quad_form(A_ext(:, i), Gamma(:, :, i)) ...
                    + (1-params.alpha)*A_ext(:, i)'*A_ext(:, i);
            end
        end

        C1 = B(A, params.d_rings, params.f_grid, params.theta_d);
        C2 = B(A, params.d_rings, params.f_grid, params.theta_BW_grid);        
        
        f0 = 1280;
        [~, idx] = min(abs(params.f_grid - f0));
        CB_bp = C2(idx:end, :);
        ref_bp = C2(end, :);
        C3 = max(max(abs(CB_bp - kron(ones(size(CB_bp, 1), 1), ref_bp)  )));
        
        minimize(obj_fun);
        subject to
            A == B_opt.'*G;
            C1 == 1;
            C2 >= params.A_BW;
            0 <= A(:) <= 1; % robustness constraint
            
%             ref_bp >= params.A_BW;
            C3 <= 0.008;
    cvx_end
%     toc
      
    costval = obj_fun;
    optval = [flip(B_opt(2:end, :)); B_opt].';
end
