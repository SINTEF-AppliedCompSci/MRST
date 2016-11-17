function [newPoints, dt] = eqInterpret(path, dt,sePtn)
% Interpolate a path with equiv distant points
% Arguments:
%   path       n*2 array of coordinates of points which are to be
%              interpolated
%   dt         distance between interpolation points
%
% Return:
%   newPoints  The interpolated points
%   dt         distance between the interpolated points

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

linesDist = sqrt(sum(diff(path,[],1).^2,2));
linesDist = [0; linesDist]; % add the starting point
cumDist = cumsum(linesDist);
s = dt*sePtn(1); e = dt*sePtn(2);
dt = (cumDist(end)-s-e)/ceil((cumDist(end)-s-e)/dt);
newPointsLoc = s:dt:cumDist(end) - e;

newPoints = interp1(cumDist, path, newPointsLoc);    
end
