function h = plotEnsembleProgress(ensemble, progress, range, h, varargin)
% Utility function for creating progress bars for running ensembles.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('width' , 800, ...
                 'height', 50 , ...
                 'onlyActive', true, ...
                 'title', '');
    opt = merge_options(opt, varargin{:});
    if nargin < 3 || isempty(h)
        figureName = 'Ensemble simulation';
        if ~isempty(opt.title)
            figureName = opt.title;
        end
        name = sprintf('%s: %s', figureName, ensemble.setup.name);
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
    title(opt.title);
end