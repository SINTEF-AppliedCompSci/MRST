function ax = drawSmoothField(ax, grid, field, samples, varargin)
    
    opt.alpha = 1;
    opt.quickclip = false;
    opt.field_threshold = [];
    opt = merge_options(opt, varargin{:});
    
    %% massaging input arguments
    if isscalar(samples)
        samples = [samples samples];
    end
    if isempty(ax)
        h = figure; plot(1);
        ax = get(h, 'currentaxes');
    end
    
    %% Establishing the continuous function
    [F, mask_fn] = continuousCelldataFunction(grid, field, opt.quickclip);
    
    %% Sampling into a regular grid
    gx = grid.cells.centroids(:,1); xvec = linspace(min(gx), max(gx), samples(1));
    gy = grid.cells.centroids(:,2); yvec = linspace(min(gy), max(gy), samples(2));
    [x, y] = meshgrid(xvec, yvec);
                      
    Fsampled = F(x, y);
    if ~opt.quickclip
        Fsampled(~mask_fn(x, y)) = NaN;
    end
    
    % also mask away any values below a specified threshold
    if ~isempty(opt.field_threshold)
        Fsampled(Fsampled < opt.field_threshold) = NaN;
    end
    
    h = pcolor(ax, xvec, yvec, reshape(Fsampled, samples(1), samples(2)));%, 'edgecolor', 'none');
    set(h, 'edgecolor', 'none');
    set(h, 'facealpha', opt.alpha);
    view(0, 90);
end
