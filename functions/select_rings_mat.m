function T = select_rings_mat(ring_idx, params)

    Mp_cumsum = cumsum(params.Nm).';
    start_idx = [1, Mp_cumsum(1:end-1)+1];
    end_idx = Mp_cumsum;

    T = zeros(params.N, params.M);
    for i = ring_idx
        T(start_idx(i):end_idx(i), i) = 1;
    end
end