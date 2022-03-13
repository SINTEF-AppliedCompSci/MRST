function map = cartBlockMap(G, rock, partition, varargin)
% Find identical coarse blocks in a partition of a Cartesian grid. This
% function returns a map, which split all coarse cells into max(map) number
% of groups. All coarse cell j with map(j)==i are in group i, and are
% considered identical, i.e. they share the exact same properties and are
% equivalent wrt. upscaling.

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

opt = struct( ...
    'satnum',    [] ...
    );
opt = merge_options(opt, varargin{:});


ncc  = numel(unique(partition)); % Number of coarse cells

% Methodology: we start with all coarse cells in the same group (group 1).
% Then we split cells into a new group if these cell do not match the other
% cells in the same group. We check for one property at a time. We will end
% up with max(map) number of groups, where all coarse cells in the same
% group share the exact same properties.

map = ones(ncc,1);

% We start by getting all fine cells for each coarse cell.
% fcells{i} contains all fine cells of coarse cell i
fcells = cellfun(@(x) find(partition == x), ...
    mat2cell((1:ncc)',ones(ncc,1)), 'UniformOutput', false);

% Group based on number of fine cells
% Get number of fine cells in each coarse cell
% nfcells = cellfun(@(x) numel(x), fcells);
% groups_nfc = unique(nfcells);
% for i = 2:numel(groups_nfc)
%    map(nfcells==groups_nfc(i)) = max(map) + 1; % Create new group
% end

% Group based on number of fine cells
nfcells = cellfun(@(x) numel(x), fcells); % # fine cells in each block
map = groupBasedOnBlockValues(map, nfcells);

% Group based on number of cell volumes
%map = groupBasedOnCellValues(map, fcells, G.cells.volumes);

% Group based on porosity
map = groupBasedOnCellValues(map, fcells, rock.poro);

% Group based on absolute permeability
for i = 1:size(rock.perm, 2)
    map = groupBasedOnCellValues(map, fcells, rock.perm(:,i));
end

% Group based on saturation regions
if ~isempty(opt.satnum)
    map = groupBasedOnCellValues(map, fcells, opt.satnum);
end


end


function map = groupBasedOnBlockValues(map, val)
% Separate new groups based on block values
% The vector val containes one value per coarse block
uniqueval = unique(val, 'stable');
for i = 2:numel(uniqueval)
    map(val==uniqueval(i)) = max(map) + 1; % Create new group
end
end


function map = groupBasedOnCellValues(map, fcells, val)
% Separate new groups based on cell values
% val is vector of size G.cells.num
ngroups = max(map);
for i = 1:ngroups
    v = cellfun(@(x) val(x), fcells(map==i), ...
        'UniformOutput', false);
    v = [v{:}]; % values of all cells in group i
    inx  = find(map==i); % cells in current group
    while true
        newgroup = [false ~all(bsxfun(@minus, v(:,2:end), ...
            v(:,1) ) == 0, 1)];
        if ~any(newgroup)
            break % we are done
        end
        v = v(:, newgroup); % Extract values of new group
        inx  = inx(newgroup);
        map(inx) = max(map) + 1; % Create new group
    end
end
end
