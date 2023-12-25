function plot_beampattern_vs_f_vs_theta(varargin)%, C

    names_idx = 1:2:length(varargin);
    values_idx = names_idx + 1;
    vargnames = [varargin{names_idx}];
    idx = find(vargnames == "bp");
    if isempty(idx)
        disp("Error plot_beampattern_vs_f_vs_theta")
    else
        bp = varargin{values_idx(idx)};
    end
    
    idx = find(vargnames == "f");
    if isempty(idx)
        disp("Error plot_beampattern_vs_f_vs_theta")
    else
        f = varargin{values_idx(idx)};
    end
    
    idx = find(vargnames == "theta");
    if isempty(idx)
        disp("Error plot_beampattern_vs_f_vs_theta")
    else
        theta = varargin{values_idx(idx)};
    end

    idx = find(vargnames == "normalize");
    if isempty(idx)
        normalize = false;
    else
        normalize = varargin{values_idx(idx)};
    end
    
    idx = find(vargnames == "plt_contour");
    if isempty(idx)
        plt_contour = false;
    else
        plt_contour = varargin{values_idx(idx)};
    end
    
    
    Ba = abs(bp);
    if normalize
        Ba = bsxfun(@rdivide, Ba, max(Ba,[],2)); % devide each row by it's max value
    end
    Ba_dB = 20*log10(Ba);

    figure;
    imagesc(theta, f, Ba_dB);
    axis 'xy'
    colorbar;
%     colormap jet;
    colormap turbo;
    caxis([-50 0]);
    xlabel('$\theta [deg]$', 'Interpreter', 'latex');
    ylabel('$ f [Hz]$', 'Interpreter', 'latex');
%     xlim([-90, 90]);

    if plt_contour
        v = [-3,-3];
        hold on;
        contour(theta, f, Ba_dB, v, 'LineColor', 'white', 'LineWidth', 2, 'LineStyle', ':')
    end
end

