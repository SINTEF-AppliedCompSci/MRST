function ax = drawContours(ax, grid, field, num_contours, varargin)

    % Most of the varargin is fed directly to the contour function, with the
    % exception of 'fill' and 'quickclip', which we extract in the lines below
    opt = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    opt.fill      = false; % just in case this option wasn't provided, set default
    opt.quickclip = true; 
    opt = merge_options(opt, varargin{:});
    do_fill   = opt.fill;
    quickclip = opt.quickclip;

    % removing fields that are not to be passed along to 'contour'
    opt = rmfield(opt, 'fill');
    opt = rmfield(opt, 'quickclip'); 
    contour_args = interleave(fieldnames(opt), struct2cell(opt));

    samples = 200;

    gx = grid.cells.centroids(:,1);
    gy = grid.cells.centroids(:,2);

    [F, mask_fn] = continuousCelldataFunction(grid, field, quickclip);
    
    [x, y] = meshgrid(linspace(min(gx), max(gx), samples), ...
                      linspace(min(gy), max(gy), samples));

    Fsampled = F(x, y);
    if ~quickclip
        % quick clipping has not been done.  Do a proper clip using the
        % returned mask (can be slow) 
        Fsampled(~mask_fn(x, y)) = NaN;  % clipping values outside grid boundary
    end
    
    if do_fill
        %colormap(plot_target, 'summer');
        fun = Fsampled/barsa;
        cvals = linspace(min(fun(:)), max(fun(:)), num_contours);
        caxis(ax, [cvals(1), cvals(end)]);
        [~, ~] = contourf(ax, x, y, fun, num_contours, contour_args{:});
        %clabel(c, ph);
        colorbar('peer', ax);
        %legend(ph, 'show');
        % colormap summer;
        % patches = findobj(ph,'-property', 'FaceAlpha');
        % set(patches, 'FaceAlpha', 0.1);
    else
        contour(ax, x, y, Fsampled, num_contours, contour_args{:});
    end
end



    