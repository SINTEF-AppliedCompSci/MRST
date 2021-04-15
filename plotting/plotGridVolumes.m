function varargout = plotGridVolumes(G, values, varargin)
% Plot partially transparent isosurfaces for a set of values
%
% SYNOPSIS:
%   plotGridVolumes(G, data)
%   plotGridVolumes(G, data, 'pn', pv,...)
%   interpolant = plotGridVolumes(...)
%
% PARAMETERS:
%   G       - Grid data structure
%
%   values  - A list of values to be plotted
%
%
% KEYWORD ARGUMENTS:
%
%   'N'     - The number of bins used to create the isosurfaces.
%
%   'min'   - Minimum value to plot. This is useful to create plots where
%             high values are visible.
%
%   'max'   - Maximum value to plot.
%
%   'mesh'  - The mesh size used to sample the interpolant. Should be a row
%             vector of length 3. Defaults to G.cartDims.
%
%   'cmap'   - Function handle to colormap function. Using different
%             colormaps for different datasets makes it possible to create
%             fairly complex visualizations.
%
%   'basealpha' - Set to a value lower than 1 to increase transparency, set
%                 it to a larger value to decrease transparency.
%
%   'binc'  - Do not call hist on dataset. Instead, use provided bins. To
%            get good results, do *not* call binc option with unique(data):
%            Ideally, binc's values should be between the unique values.
%
%   'patchn' - Maximum number of patch faces in total for one call of
%            plotGridVolumes. If this number is large, the process may be
%            computationally intensive.
%
%   'interpolant' - If you are plotting the same dataset many times, the
%             interpolant can be returned and stored.
%
%   'extrudefaces' - Let the cell values be extrapolated to the edges of the
%                  domain. Turn this off if you get strange results.
%
% RETURNS:
%   interpolant - See keyword argument of same name.
%
% SEE ALSO:
%   `plotCellData`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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


   opt = struct('N',    20, ...
                'min', [], ...
                'max', [], ...
                'mesh', [], ...
                'cmap',  @jet,...
                'basealpha', 1, ...
                'binc',    [],...
                'patchn', 1e6,...
                'interpolant', [], ...
                'extrudefaces', true);

    opt = merge_options(opt, varargin{:});
    N = opt.N;

    assert(G.griddim == 3, 'plotGridVolumes is only supported in 3D!');

    if isempty(opt.mesh)
        mesh = G.cartDims + 2*opt.extrudefaces;
    else
        mesh = opt.mesh;
    end

    v = values;
    if ~isempty(opt.min)
        v = v(v >= opt.min);
    end

    if ~isempty(opt.max)
        v = v(v <= opt.max);
    end

    if isempty(opt.binc)
        [binc, binc] = hist(v, N);                                    %#ok
    else
       binc = opt.binc;
       N = numel(binc);
    end


    gc = G.cells.centroids;
    bf = boundaryFaces(G);
    gfc = G.faces.centroids(bf, :);
    if opt.extrudefaces
        cellfa = sum(G.faces.neighbors(bf,:),2);
        gc = [gc;  gfc];
        values = [values; values(cellfa)];
    end

    Mgc = max(gc); Mgc = min([Mgc; max(G.nodes.coords)]);
    mgc = min(gc); mgc = max([mgc; min(G.nodes.coords)]);

    tmp = cell(G.griddim, 1);
    for i = 1:G.griddim
        tmp{i} = linspace(mgc(i),Mgc(i),mesh(i)+1);
    end

    [x,y,z] = meshgrid(tmp{:});

    if (ispc || ismac) && ~strcmpi(get(gcf, 'Renderer'), 'OpenGL')
        set(gcf, 'Renderer', 'OpenGL');
    end

    c = opt.cmap(N);
    if isempty(opt.interpolant)
        if exist('scatteredInterpolant', 'file')
            F = scatteredInterpolant(gc(:,1), gc(:,2), gc(:,3), values);
        elseif exist('TriScatteredInterp', 'file')
            F = TriScatteredInterp(gc(:,1), gc(:,2), gc(:,3), values); %#ok<DTRIINT>
        else
            assert (exist('griddata3', 'file') == 2, ...
                    'Function GRIDDATA3 does not exist?');
            F = @(xi, yi, zi) ...
                griddata3(gc(:,1), gc(:,2), gc(:,3), values, xi, yi, zi);
        end
        interpolant = F(x,y,z);
    else
        interpolant = opt.interpolant;
    end
    for i = 1:N
        p = patch(isosurface(x, y, z, interpolant, binc(i)), ...
                  'FaceColor', c(i,:), 'EdgeColor', 'none', ...
                  'FaceAlpha', opt.basealpha*i/(2*N));
        fn = size(get(p, 'Faces'),1);
        set(get(p, 'Parent'), 'ZDir', 'reverse')
        if fn > opt.patchn/N
            reducepatch(p, ceil(opt.patchn/N));
        end
    end

    if nargout
        varargout{1} = interpolant;
    end
end
