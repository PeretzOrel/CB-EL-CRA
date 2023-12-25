function plot_array_geometry(varargin)%, C

    names_idx = 1:2:length(varargin);
    values_idx = names_idx + 1;
    vargnames = [varargin{names_idx}];
    idx = find(vargnames == "elem_pos");
    if isempty(idx)
        disp("Error plot_array_geometry2")
    else
        elem_pos = varargin{values_idx(idx)};
    end
    
    idx = find(vargnames == "weights");
    if isempty(idx)
        color_by_weights = false;
    else
        color_by_weights = true;
        weights = varargin{values_idx(idx)};
    end
    
    idx = find(vargnames == "partition");
    if isempty(idx)
        partition = 1;
        color_by_partition = false;
    else
        color_by_partition = true;
        partition = varargin{values_idx(idx)};
    end

    figure;
    sp=axes;
    
    X = elem_pos(:, 1);
    Y = elem_pos(:, 2);
    Z = elem_pos(:, 3);
    
    if color_by_weights
        s = scatter3(sp, X, Y, Z, 20, weights, 'o', 'filled');
        colormap(jet);
        colorbar;
    elseif color_by_partition
        active_elem = sum(partition, 2) > 0;
        
        ii = 1;
        for i = 1:size(partition, 2)
            if sum(partition(:, i)) > 0
                partition(:, i) = partition(:, i)*ii;
                ii = ii + 1;
            end
        end
        partition = sum(partition, 2);
        s = scatter3(sp, X(active_elem), Y(active_elem), Z(active_elem), 20, partition(active_elem), 'o', 'filled');
        colormap(jet);
%         colormap(hsv);
%         colormap(parula);


        hold on;
        s = scatter3(sp, X(~active_elem), Y(~active_elem), Z(~active_elem),20,'o', 'filled','MarkerFaceColor',	[1 1 1]);    

    else
        s = scatter3(sp, X, Y, Z,20,'o', 'filled','MarkerFaceColor','k');       
    end
    
    axHdl = get(s,'parent');
    set(axHdl,'DataAspectRatio',[1 1 1]);
    view_az = 135; view_el = 20;
    view([view_az view_el]);
    axis vis3d
    axis off;
    hold on;
    zHdl = zoom;
    zHdl.setAxes3DPanAndZoomStyle(axHdl,'camera');
    
    % Get current positions of axes
    xmax = max(max(abs(X)));
    ymax = max(max(abs(Y)));
    zmax = max(max(abs(Z)));
    maxdata = max([xmax ymax zmax]);
    xmax = max(maxdata/4,xmax);  % ensure the axis span is at least maxdata/4
    ymax = max(maxdata/4,ymax);
    zmax = max(maxdata/4,zmax);
    if isfinite(xmax) && xmax~=0
        set(get(s,'parent'),'XLim',1.3*[-xmax xmax]);
    end
    if isfinite(ymax) && ymax~=0
        set(get(s,'parent'),'YLim',1.3*[-ymax ymax]);
    end
    if isfinite(zmax) && zmax~=0
        set(get(s,'parent'),'ZLim',1.3*[-zmax zmax]);
    end

    axisfactor = 1.2;
    XPos=axisfactor*xmax;
    YPos=axisfactor*ymax;
    ZPos=axisfactor*zmax;
    
    % Create pseudo axes and mark ticks
    plot3( [0,XPos],[0,0],[0,0],'k','LineWidth',1.5,'Tag','XAxis' );
    text(1.3*XPos,0,0,["$x$", "$\varphi=0$", "$\theta=90$"],'Tag','Az0El0Label','Interpreter', 'latex');
    plot3( [0,0],[0,YPos],[0,0],'k','LineWidth',1.5,'Tag','YAxis' );
    text(0,1.15*YPos,0,["$y$", "$\varphi=90$", "$\theta=90$"],'Tag','Az90El0Label','Interpreter', 'latex');
    plot3( [0,0],[0,0],[0,ZPos],'k','LineWidth',1.5,'Tag','ZAxis' );
    text(0,0,1.15*ZPos,["$z$", "$\varphi=0$", "$\theta=0$"],'Tag','Az0El90Label','Interpreter', 'latex');

    hold off;
    
    view(0,90)
end