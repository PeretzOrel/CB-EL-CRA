function plot_beampattern_vs_theta(varargin)

    names_idx = 1:2:length(varargin);
    values_idx = names_idx + 1;
    vargnames = [varargin{names_idx}];
    idx = find(vargnames == "hax");
    if isempty(idx)
        figure; 
        hax=axes;
    else
        hax = varargin{values_idx(idx)};
    end
    
    idx = find(vargnames == "bp");
    if isempty(idx)
    else
        B = varargin{values_idx(idx)};
    end
    
    idx = find(vargnames == "theta");
    if isempty(idx)
        theta = linspace(-90, 90, 2000);
    else
        theta = varargin{values_idx(idx)};
    end

    
    Ba = abs(B);
%     Ba = bsxfun(@rdivide, Ba, max(Ba,[],2)); % devide each row by it's max value
    Ba_dB = 20*log10(Ba);
    plot(hax, theta, Ba_dB);
    xlabel(hax, '$\theta [deg]$', 'Interpreter', 'latex');
    ylabel(hax, '$ Mag [dB]$', 'Interpreter', 'latex');
    ylim(hax, [-60 2]);
    xlim(hax, [min(theta) max(theta)]);
    grid(hax, 'on');
    hold(hax, 'on');
end

