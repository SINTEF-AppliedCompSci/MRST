function p = translateLine(points,a)
% translateLine(points,a) translates the line given by an n-by-2 vector of
% coincident 'points' by a normal distance 'a' in 2D space to create a
% rectangle given by the 2*n-by-2 vector 'p'.

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


endp = [points(1,:);points(end,:)];
diff1 = diff(endp,1);
m = diff1(2)/diff1(1); % slope
c = endp(1,2)-m*(endp(1,1)); % c = y1 - m*x1
xt = zeros(size(points,1),1); yt = xt;
for i = 1:size(points,1)
    px = points(i,1); py = points(i,2);
    xt(i) = (py + px/m - a*sqrt(1+m^2) - c)/(m+1/m);
    yt(i) = m*xt(i) +  a*sqrt(1+m^2) + c;
end
p = [points;xt,yt];
return