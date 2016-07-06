function dist = pdist_euclid(x)
% Pairwise euclidian distance between pairs of objects.

assert(~isempty(x) && numel(size(x))<=2, 'x must be a 1 or 2 dimensional non-empty matrix');
ind = nchoosek(1:size(x,1),2);
Xi = ind(:,1);
Yi = ind(:,2);
X = x';
diff = X(:,Xi) - X(:,Yi);
if numel(diff) == 2
    dist = norm(diff);
else
    dist = sqrt(sumsq(diff,1));
end
return

function s = sumsq(x,dim) % 2D matrices only
if dim == 1
    s = zeros(1,size(x,2));
    for i = 1:size(x,2)
        s(i) = norm(x(:,i))^2;
    end
elseif dim == 2
    s = zeros(size(x,1),1);
    for i = 1:size(x,1)
        s(i) = norm(x(i,:))^2;
    end
else
    error('Input must be a 2D matrix and dim must be 1 or 2');
end
return

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

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