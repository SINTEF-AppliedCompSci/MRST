function CG = setupMexInteractionMapping(CG, varargin)
% Setup mex mapping for cppMultiscaleBasis

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

    opt = struct('verbose', mrstVerbose());
    opt = merge_options(opt, varargin{:});
    
    if isfield(CG.cells, 'support_mex')
        return
    end
    
    CG.cells.support_mex = getGridData(CG, opt);
end

function mexSupport = getGridData(CG, opt)
    G = CG.parent;

    interaction = CG.cells.interaction;
    bnd = cell(CG.cells.num, 1);

    isGlobalBound = false(G.cells.num, 1);
    ticif(opt.verbose);
    for i = 1:CG.cells.num
        dispif(opt.verbose, 'Processing block %d of %d\n', i, CG.cells.num);
        c = interaction{i};
        % Get the boundary cells - cells that will never have value > 0.
        [slf, othr]= boundaryNeighbors(G, c);
        bnd{i} = othr;
        isGlobalBound(othr) = true;
    end
    tocif(opt.verbose);
    % Tag global boundary as 2
    fn1 = @(x) cellfun(@(y) 2*isGlobalBound(y), x, 'UniformOutput', false);
    % Tag local boundary as 1
    fn2 = @(x) cellfun(@(y) 0*y + 1, x, 'UniformOutput', false);

    cells = horzcat(interaction, bnd)';
    indicator = horzcat(fn1(interaction), fn2(bnd))';

    [sort_cells, sort_ind] = deal(cell(numel(interaction), 1));
    ticif(opt.verbose);
    for i = 1:CG.cells.num
        dispif(opt.verbose, 'Sorting block %d of %d\n', i, CG.cells.num);
        [sort_cells{i}, ix] = sort(vertcat(cells{:, i}));
        tmp = vertcat(indicator{:, i});
        sort_ind{i} = tmp(ix);
    end
    tocif(opt.verbose);
    
    cm = [vertcat(sort_cells{:}), vertcat(sort_ind{:})];
    counts = sum(cellfun(@numel, cells), 1);
    offsets = cumsum([1; counts(:)]);

    I = zeros(size(cm, 1), 1);
    for i = 1:CG.cells.num
        subs = offsets(i):(offsets(i+1)-1);
        I(subs) = CG.partition(cm(subs, 1)) == i;
    end

    % Convert to int
    int_offsets = toIntegerIndex(offsets);
    support = toIntegerIndex(cm(:, 1));
    types = toIntegerValue(cm(:, 2));
    
    mexSupport = struct('offsets', int_offsets,...
                        'celltypes', types, ...
                        'support', support, ...
                        'sorted_cell', {sort_cells});
end

function v = toIntegerIndex(v)
    v = toIntegerValue(v - 1);
end

function v = toIntegerValue(v)
    v = int32(v);
end