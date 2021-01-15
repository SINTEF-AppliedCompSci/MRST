function m = flood(G, start, faceStatus)
%
%
% SYNOPSIS:
%   function m = flood(G, start, faceStatus)
%
% DESCRIPTION:
%
% PARAMETERS:
%   G           - Grid
%   start       - Coordinate for starting point
%   faceStatus  - Face status vector from sliceGrid. If calling
%                 sliceGrid as [G, gix] = sliceGrid(...), then
%                 faceStatus = gix.new.faces.
%
% RETURNS:
%   m - Vector with flooded cells marked 1
%
% EXAMPLE:
%
% SEE ALSO: sliceGrid.m
%

    m = zeros(G.cells.num, 1);
    checked = zeros(G.cells.num, 1);
    notDone = findEnclosingCell(G, start);
    cnt = 0;
    maxCnt = intmax;

    % The sliceGrid source fixes the value of the faces forming the
    % cutting plane
    faces_along_plane = 3;

    while ~isempty(notDone) & cnt < maxCnt
        % Pop first element
        e = notDone(1);
        notDone(1) = [];
        m(e) = 1;
        checked(e) = 1;

        % Exclude faces that form the cut plane
        f = G.cells.faces(G.cells.facePos(e):G.cells.facePos(e+1)-1, :);
        f(faceStatus(f) == faces_along_plane) = [];

        % Add neighboring cells, excluding 0, e and checked elements
        ne = G.faces.neighbors(f, :);
        ne = ne(:);
        ne(ne == 0) = [];
        ne(ne == e) = [];
        ne(checked(ne) == 1) = [];
        notDone = unique([notDone; ne]);

        % Counter
        cnt = cnt+1;
    end

    if cnt == maxCnt
        error('flood failed');
    end

end
