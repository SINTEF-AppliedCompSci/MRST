function m = markCutGrids(G, faceStatus, varargin)
%
%
% SYNOPSIS:
%   function m = markCutGrids(G, faceStatus)
%
% DESCRIPTION:
%
% PARAMETERS:
%   G               - Grid
%   faceStatus      - Face status vector from sliceGrid. If calling
%                     sliceGrid as [G, gix] = sliceGrid(...), then
%                     faceStatus = gix.new.faces.
%   facesAlongPlane - Optional number indicating faceStatus value of
%                     the internal boundary. Default is 3.
%   legacy          - Optional force of older code without graph.
%   start           - Legacy code requires starting position
%
% RETURNS:
%   m - Vector with marked cells
%
% EXAMPLE:
%
% SEE ALSO: sliceGrid.m
%

    % The sliceGrid source fixes the value of the faces forming the
    % cutting plane
    opt = struct('facesAlongPlane', 3, ...
                 'legacy', false, ...
                 'start', []);
    opt = merge_options(opt, varargin{:});

    if verLessThan('matlab', '8.6') | opt.legacy
        % Flood
        warning('This is slow')
        m = zeros(G.cells.num, 1);
        checked = zeros(G.cells.num, 1);
        notDone = findEnclosingCell(G, opt.start);
        cnt = 0;
        maxCnt = intmax;

        while ~isempty(notDone) & cnt < maxCnt
            % Pop first element
            e = notDone(1);
            notDone(1) = [];
            m(e) = 1;
            checked(e) = 1;

            % Exclude faces that form the cut plane
            f = G.cells.faces(G.cells.facePos(e):G.cells.facePos(e+1)-1, :);
            f(faceStatus(f) == opt.facesAlongPlane) = [];

            % Add neighboring cells, excluding 0, e and checked elements
            ne = G.faces.neighbors(f,:);
            ne = ne(:);
            ne(ne == 0) = [];
            ne(ne == e) = [];
            ne(checked(ne) == 1) = [];
            notDone = [notDone; ne];

            % Counter
            cnt = cnt+1;
        end

        if cnt == maxCnt
            error('markCutGrids failed');
        end

    else
        % Use graph
        faces = (1:G.faces.num)';
        faces(faceStatus == opt.facesAlongPlane) = [];

        % Get neighboring cells
        neighs = G.faces.neighbors(faces, :);
        idx = all(neighs > 0, 2);
        neighs = neighs(idx, :);

        % Form adjacency matrix
        g = graph(neighs(:,1), neighs(:,2));
        m = conncomp(g)';

    end

end
