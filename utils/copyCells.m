function [G, cno] = copyCells(G, cells, varargin)
%Copy a set of cells in a grid
%
% SYNOPSIS:
%   [G, cno] = coypCells(G, cells);
%   [G, cno] = coypCells(G, cells, 'pn1', pv1, ...);
%
% DESCRIPTION:
%   This function copies a set of cells and adds them to the grid. The
%   resulting grid will have G.cells.num + numel(cells) cells, and the
%   copied cells with share the same faces as their originals. Connections
%   are not copied, and connections from the new cells to any of the cells
%   in the new grid can be specified as non-neighboring connections.

% REQUIRED PARAMETERS:
%   G     - Grid structure as described by grid_structure.
%   cells - Cells to be copied, either as a vector of cell numbers, or a
%           logical mask onto the original grid cells.
%
% OPTIONAL PARAMETERS:
%   nnc - Non-neighboring connections for the new cells. A struct
%         with a field `cells` with an mx2 matrix of connections for copied
%         cells (and optionally `trans`), or an mx2 matrix of connections
%         for copied cells. Values in the first column of of the connection
%         matrix are interpreted as new cell indices, given by the cell
%         number of its original, and values in the second column are
%         interpreted as cells in the original grid. If `trans` is not
%         provided, the corresponding rows in G.nnc.trans will be populated
%         with nans. Default value: [].
%
% RETURNS:
%   G   - Grid with copied cells included (with nncs if given)
%   cno - Vector of cell numbers for the new cells
%
% EXAMPLES:
%   Copy cell number 1 of a grid and add an nnc between the copy and the
%   original cell:
%
%   G = cartGrid([2,2]);
%   G = copyCells(G, 1, 'nnc', [1,1]);
%

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    % Process input
    %---------------------------------------------------------------------%
    % Optional arguments
    opt = struct('nnc', []); % Connection cells
    opt = merge_options(opt, varargin{:});
    % Validate cells to be copied
    [cells, opt] = validateInput(G, cells, opt);
    %---------------------------------------------------------------------%
    
    % Copy cells
    %---------------------------------------------------------------------%
    % Update number of cells and compute cellnumbers for new cells
    n0          = G.cells.num;
    G.cells.num = n0 + numel(cells);
    cno         = (1:numel(cells))' + n0;
    
    % Update facePos for new cells
    nf              = diff(G.cells.facePos); nf = nf(cells);
    facePos         = cumsum([G.cells.facePos(end); nf]);
    G.cells.facePos = [G.cells.facePos; facePos(2:end)];
    % Set faces for new cells
    faces = G.cells.faces(mcolon(G.cells.facePos(cells), ...
                                 G.cells.facePos(cells + 1) - 1),:);
    G.cells.faces = [G.cells.faces; faces];
    
    if any(strcmpi(G.type, 'computegeometry'))
        % Copy geometry for new cells
        G.cells.centroids = G.cells.centroids([(1:end)'; cells],:);
    end
    %---------------------------------------------------------------------%

    % Add connections
    %---------------------------------------------------------------------%
    if ~isempty(opt.nnc)
        % Add connections for new cells
        if ~isfield(G, 'nnc'), G.nnc = struct('cells', [], 'trans', []); end
        % Construct map from oiriginals to copies
        cmap = nan(G.cells.num, 1); cmap(cells) = cno;
        G.nnc.cells = [G.nnc.cells; [cmap(opt.nnc.cells(:,1)), opt.nnc.cells(:,2)]];
        G.nnc.trans = [G.nnc.trans; opt.nnc.trans];
        assert(~any(isnan(G.nnc.cells(:))), ['Invalid connection ', ...
            'matrix for nncs! See documentation for correct use.']);
    end
    %---------------------------------------------------------------------%
    
end

%-------------------------------------------------------------------------%
function [cells, opt] = validateInput(G, cells, opt)

    % Check cells to be copied
    if islogical(cells), cells = find(cells); end;
    cells = reshape(cells, [], 1);
    assert(min(cells) >= 1 && max(cells) <= G.cells.num,             ...
        ['Input cells must correspond to cell numbers in the range', ...
         '[1, G.cells.num]']                                       );
    
    if isempty(opt.nnc), return; end
    % Check nncs
    if isnumeric(opt.nnc), opt.nnc = struct('cells', opt.nnc); end
    if ~isfield(opt.nnc, 'trans'), opt.nnc.trans = nan(size(opt.nnc.cells,1), 1); end
    
end