function pick = findClosestCell(g,x,y,sub)
% Find the cell with center closest to the point (x,y)
%
% SYNOPSIS
% pick = findClosestCell(g,x,y,sub)
%
% PARAMETERS
%       g - mrst grid structure
%       x,y - x and y coordinates
%
%       OPTIONAL
%       sub - a subset of the cells to search within
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.

% search among all cells if no restriction are given
if nargin <4
    sub = 1:g.cells.num;
end

% pick out the centroids
xc = g.cells.centroids(sub,1);
yc = g.cells.centroids(sub,2);

% find the closest cell, one point at the time
pick = zeros(numel(x),1);
for i = 1:numel(x)
    d = sqrt((xc - x(i)).^2 + (yc - y(i)).^2);

    [~,ind] = min(d);
    if nargin == 4
        if numel(sub)==g.cells.num
            sub = find(sub);
        end
        ind = sub(ind);
    end
    pick(i) = ind;
end