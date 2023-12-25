function plot_weights(H, params)
    figure;
    imagesc(1:params.M, params.f_grid, abs(H.'));
    axis xy; 
    colorbar; 
    colormap turbo;
    ax = gca;
    new_label = arrayfun(@(x) num2str(x/1000), ax.YTick, 'UniformOutput', false);
    ax.YTickLabel = new_label;
    % new_label = arrayfun(@(x) num2str(x), 1:M, 'UniformOutput', false);
    % ax.XTickLabel = new_label;
    % caxis([0, 1]);

    figure;
    for i = 1:params.M
        plot(params.f_grid, abs(H(i, :)),'DisplayName',"Ring #: " + string(i)); hold on;
        legend
        grid on;
    end

end