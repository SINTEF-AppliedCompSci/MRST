function Gr = readAtlasGrids(varargin)

if nargin > 0
    gdir = varargin{1};
else
    gdir = fullfile(VEROOTDIR, 'data', 'atlas');
end
dir_grid = dir(gdir);
grids = {dir_grid(~[dir_grid.isdir]).name};

Gr = {};
for i = 1:numel(grids)
    g = grids{i};
    % Skip projections
    if strcmpi(g(end-3:end), '.prj')
        continue
    end
    [meta, data] = readAAIGrid(fullfile(gdir, g));

    G = processAAIGrid(meta, data, true, []);
    tmp = regexp(g, '_', 'split');
    G.name = tmp{1};

    % Designate the type based on the file names
    ind = strcmpi(tmp, 'thickness') | strcmpi(tmp, 'top');
    if any(strcmpi(tmp, 'thickness'))
        G.variant = 'thickness';
    else
        G.variant = 'top';
    end
    G.name = [tmp{~ind}];
    G.meta = meta;
    G.data = data;
    Gr{end+1} = G;

    assert(isfield(G.cells,'cellNodes'))
end
end
