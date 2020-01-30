function [p, removed] = surfaceSufCond2D(p, F)
% Enforces the sufficient surface condition.
%
% SYNOPSIS:
%   [p, removed] = surfaceSufCond2D(p,F)
%
% PARAMETERS:
%   p         - A n X 2 array of points
%   F         - A surface struct as returned from surfaceSites2D(..)
%
% RETURNS: 
%   p         - All sites from the input p that does not violate the
%               sufficient surface condition, that is, all sites from input
%               p that lie innside the circles in F are removed.
%   removed   - A logical array of length n. If element removed(i) is true
%               the site p(i,:) from input p violates the sufficient
%               surface condition.
% 
% EXAMPLE:
%   p = rand(200,2);
%   F = surfaceSites2D({[0.2,0.5;0.8,0.5]},0.2);
%   [~,r] = surfaceSufCond2D(p,F);
%   theta = linspace(0,2*pi)';
%   figure(); hold on
%   for i = 1:size(F.c.CC,1)
%     X = repmat(F.c.CC(i,:),100,1) + F.c.R(i)*[cos(theta), sin(theta)];
%     plot(X(:,1),X(:,2),'k');
%   end
%   plot(p(r,1), p(r,2),'.r');
%   plot(p(~r,1), p(~r,2),'.b');
%
% SEE ALSO
%   compositePebiGrid2D, pebi, pebiGrid, clippedPebi2D, surfaceSufCond3D

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2015-2020 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

TOL = 10 * (max(F.f.pts(:)) - min(F.f.pts(:))) * eps;
nc = size(F.c.CC,1);
np = size(p,1);

CRSqr = F.c.R.^2;
removed = zeros(np,1);
for i = 1:nc
  distSqr = sum(bsxfun(@minus, F.c.CC(i,:), p).^2,2);
  removed = removed + (distSqr<CRSqr(i)-TOL);
end
removed = logical(removed);
p = p(~removed,:);
end
