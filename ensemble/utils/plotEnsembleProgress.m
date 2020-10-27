function h = plotEnsembleProgress(ensemble, progress, range, h, varargin)
    opt = struct('width' , 800, ...
                 'height', 50 , ...
                 'onlyActive', true);
    opt = merge_options(opt, varargin{:});
    if nargin < 3 || isempty(h)
        name = sprintf('Ensemble simulation: %s', ensemble.setup.name);
        h = figure('Name', name);
    end
    clf(h);
    
    total = numel(progress);
    done = isinf(progress);
    active = progress > 0 & progress < inf;
    keep = active;
    range = range(keep);
    progress = progress(keep);
    n = numel(range)+1;
    height = 0.9/n;
    colors = lines(2);
    for i = 1:n-1
        str = sprintf('Ensemble member %d', range(i));
        simpleUIbar(h, min(progress(i), 1), (i-1)*height, height, str, 'FaceColor', colors(1,:));
    end
    str = sprintf('Total progress (%u/%u finished)', nnz(done), total);
    simpleUIbar(h, (nnz(done) + sum(progress))/total, (n-1)*height, height, str, 'FaceColor', colors(2,:));
    
    h.Position(3) = opt.width;
    h.Position(4) = opt.height*n;
end