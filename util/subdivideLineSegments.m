function [Pts, origIx, ptsUsed] = subdivideLineSegments(path, dt, sePtn)
% Subdivide line segments with almost equidistant points, but respect input
% Arguments:
%   path       n*2 array of coordinates defining the line segments
%
%   dt         default distance for subdivision points
%
%   sePtn      sePtn*dt defines increments used to trim the curve at the
%              start and end points 
%
% Return:
%   Pts        the new points
%
%   origIx     which among the points Pts are points from input path 
%
%   ptsUsed    index of the end points of the line segments we use to
%              define the sample point (if sePtn is positive, we may have
%              trimmed away line segments at each end)

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2019 Knut-Andreas Lie. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

sDist = sqrt(sum(diff(path,[],1).^2,2));
cCoor = [0; cumsum(sDist)];
split = ceil(sDist./dt);
nCoor = nan(sum(split)+1,1);

sPt  = sePtn(1)*dt;
ePt  = cCoor(end)-sePtn(2)*dt;
n    = 1;
flag = false(numel(sDist)+1,1);
for i=1:numel(sDist)
    if (cCoor(i+1)<sPt) || (cCoor(i)>ePt)
        split(i)=0; continue
    end
    ds          = min(cCoor(i+1),ePt)-max(cCoor(i),sPt);
    split(i)    = ceil(ds/dt);
    nCoor(n:n+split(i)) = max(cCoor(i),sPt):ds/split(i):min(cCoor(i+1),ePt);
    n           = n + split(i);
    flag(i:i+1) = true;
end
nCoor   = nCoor(~isnan(nCoor));
split   = split(split>0);
origIx  = [1; cumsum(split)+1];
origIx  = origIx(1+(sePtn(1)>0):end-(sePtn(2)>0));
Pts     = interp1(cCoor, path, nCoor);
ptsUsed = find(flag); 
end