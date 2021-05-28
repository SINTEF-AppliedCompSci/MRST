function model = addCoarseOperatorsMultiVE(model)
% Add various extra coarse-scale information for a multi-VE model
%
% SYNOPSIS:
%   model = addCoarseOperatorsMultiVE(model)
%
% REQUIRED PARAMETERS:
%   model - Output from generateCoarseGridMultiVE
%
% OPTIONAL PARAMETERS:
%   None.
%
% RETURNS:
%   model - A modified model.
%
% SEE ALSO:
%   addCoarseOperatorsMultiVE

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

    G = model.G;
    if ~isfield(G, 'partition') || ~isfield(G.cells, 'discretization')
        error('Function only applicable to coarse hybrid VE grids from generateCoarseGridMultiVE');
    end
    N = model.operators.N;
    n1 = N(:, 1);
    n2 = N(:, 2);

    d = G.cells.discretization;
    % Change in discretization type
    transition = d(n1) ~= d(n2);

    % Internal connections inside same VE zone
    model.operators.connections.veInternalConn = ...
        d(n1) > 1 & ~transition;
    % Between two fine cells
    model.operators.connections.fineInternalConn =...
        d(n1) == 1 & ~transition;
    % Between two VE zones horizontally
    ve_trans = transition & d(n1) > 0 & d(n2) > 0;
    diff_column = G.cells.columns(n1) ~= G.cells.columns(n2);
    model.operators.connections.veTransitionHorizontalConn = ...
        ve_trans & diff_column;
    % Between two VE zones vertically
    model.operators.connections.veTransitionVerticalConn = ...
        ve_trans & ~diff_column;
    % From one VE zone to fine
    model.operators.connections.veToFineConn = ...
        transition & (d(n1) == 1 | d(n2) == 1);

    globalCoarseFace = find(all(G.faces.neighbors > 0, 2));

    nt = numel(globalCoarseFace);

    transFaceTop = zeros(nt, 2);
    transFaceBottom = zeros(nt, 2);
    faceHeights = zeros(nt, 2);
    for i = 1:nt
        f = globalCoarseFace(i);
        fine = G.faces.fconn(G.faces.connPos(f):G.faces.connPos(f+1)-1);
        C = G.faces.neighbors(f, :);

        c0 = G.parent.faces.neighbors(fine, :);
        blockNo = G.partition(c0);
        if numel(fine) == 1
            blockNo = blockNo';
        end
        c = c0;

        sorted = blockNo(:, 1) == C(1);
        c(:, 1) = sorted.*c0(:, 1) + ~sorted.*c0(:, 2);
        c(:, 2) = sorted.*c0(:, 2) + ~sorted.*c0(:, 1);
        if numel(fine) == 1
            t = G.parent.cells.topDepth(c)';
            b = G.parent.cells.bottomDepth(c)';
        else
            t = min(G.parent.cells.topDepth(c));
            b = max(G.parent.cells.bottomDepth(c));
        end
        H = G.cells.height(C)';
        B = G.cells.bottomDepth(C)';
        T = G.cells.topDepth(C)';

        if 1
            % Think this is ok - scale the height by the fraction.
            % Sort of assumes that columns are vector-like with no
            % deviations.
            faceHeights(i, :) = H.*(b - t)./(B - T);
        else
            % "Real geometry" version - but less correct for the VE maths
            faceHeights(i, :) = b - t;
        end
        transFaceTop(i, :) = t;
        transFaceBottom(i, :) = b;
    end
    % Top of the columns connected to the face
    model.operators.connections.faceTopDepth = transFaceTop;
    % Bottom of the columns connected to the face
    model.operators.connections.faceBottomDepth = transFaceBottom;
    % Height of the cell-column connected to the specific face
    model.operators.connections.faceHeight = faceHeights;
end