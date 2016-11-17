function [p, removed] = faultSufCond(p, F)
% Enforces the sufficient fault condition.
%
% SYNOPSIS:
%   [p, removed] = faultSufCond(p,F)
%
% PARAMETERS:
%   p         - A n X 2 array of points
%   F         - A fault struct as returned from createFaultGridPts(..)
%
% RETURNS: 
%   p         - All points from the input p that does not violate the
%               sufficient fault condition, that is, all points from input
%               p that lie innside the circles in F are removed.
%   removed   - A logical array of length n. If element removed(i) is true
%               the point p(i,:) from input p violates the sufficient fault
%               condition.
% 
% EXAMPLE:
%   p = rand(200,2);
%   F = createFaultGridPoints({[0.2,0.5;0.8,0.5]},0.2);
%   [~,r] = faultSufCond(p,F);
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
%   compositePebiGrid, pebi, pebiGrid, clippedPebi2D, faultSufCond3D

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

TOL = 50*eps;
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