function [p,offSet] = ballInt(CC1, CR1, CC2, CR2, CC3, CR3)
% Calculates the intersection of three balls. 
% 
% Supports intersection of mulitiple sets of balls. 
%
% SYNOPSIS:
%   [p, offSet = ballInt(CC1, CR1, CC2, CR2, CC3, CR3)
%
% PARAMETERS:
%   CC1       - A nx3 array containing the ball center of balls number 1.
%   CR1       - A nx1 array containing the radii of balls number 1.
%   CC2       - A nx3 array containing the ball center of balls number 2.
%   CR2       - A nx1 array containing the radii of balls number 2.
%   CC3       - A nx3 array containing the ball center of balls number 3.
%   CR3       - A nx1 array containing the radii of balls number 3.
%
% RETURNS:
%   p         - A 2nx3 array containing the intersection of the balls.
%               p(2*i-1:2*i,:) contains the intersection of balls CC1(i,:),
%               CC2(i,:) and CC3(i,:).
%   offSet    - A 2nx1 array containing the distance between intersections.
%               If i is odd, offset(i) is the distance between p(i,:) and
%               p(i+1,:). If i is even, offset(i) is the distance between
%               p(i,:) and p(i-1,:).
%
% EXAMPLE:
% CC1 = rand(5,3);
% CC2 = rand(5,3);
% CC3 = rand(5,3);
% CR  = ones(5,1);
% 
% p = ballInt(CC1, CR, CC2, CR, CC3, CR)
%
% SEE ALSO:
%   surfaceSites3D, calcNormals, isColinear

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

if isempty(CC1) || isempty(CC2)
  p = [];
  return
end

assert(all(size(CC1)==size(CC2) & size(CC1)==size(CC3)),...
  'The number of circle centers must be the same for all circle sets');
assert(size(CC1,1)==size(CR1,1) ...
    && size(CC2,1)==size(CR2,1) ...
    && size(CC3,1)==size(CR3,1),...
    'The number of circle centers must equal the number of radii');

dim = size(CC1,2);

% Create basis vector pointing from c1 to c2.
p1 = CC2 - CC1;
dx = sqrt(sum((p1).^2,2));               % Distance between c1 & c2
nx = bsxfun(@rdivide,p1,dx);             % basis vector 1

% Create the basis vector - in the plane defined by c1, c2 and c3 - that is
% orthogonal to nx.
p2 = CC3 - CC1;
midX = dot(nx, p2,2);
ny   = p2 - bsxfun(@times,midX,nx);
ny   = bsxfun(@rdivide, ny, sqrt(sum(ny.^2,2)));

% Create basis vector orthogonal to nx and ny.
nz = cross(nx,ny);

% Find coordinates in new basis
X = (dx.^2 - CR2.^2 + CR1.^2)./(2*dx);
midY = dot(ny,p2,2);
Y = (CR1.^2 - CR3.^2 + midY.^2 + midX.^2)./(2*midY) - midX./midY.*X;
Z  = sqrt(CR1.^2 - X.^2-Y.^2);

% Transform back to cartesian coordinates.
top = CC1 +   bsxfun(@times, X, nx)  ...
          +   bsxfun(@times, Y, ny)  ...
          +   bsxfun(@times, Z, nz);
bot = top - 2*bsxfun(@times, Z, nz);

% Put together result
p      = reshape([top,bot]',dim,[])';
offSet = reshape([Z,Z]',1,[])';
end
