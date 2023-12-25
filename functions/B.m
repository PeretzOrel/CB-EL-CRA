function bp = B(H, func, f, theta)
    [theta_mat, f_mat] = meshgrid(theta, f);
    d_f_theta = cell2mat(arrayfun(func, f_mat(:), theta_mat(:) ,'UniformOutput', false).');
    d_f_theta = reshape(d_f_theta, [], length(f), length(theta));
    HH = repmat(H, 1, 1, length(theta));
    bp = HH.*d_f_theta;
    bp = squeeze(sum(bp, 1));
end
