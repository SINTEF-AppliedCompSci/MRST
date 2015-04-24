doprint = true;

%% Example of MaxMin Criterion
clf;
x = [.2 1 1.5 .7 -0.2 NaN];
y = [0 .1 0.5 1    .8 NaN];
subplot(1,2,1);
plot(x,y,'o','MarkerSize',8);
I = [1 2 5 1 6 2 3 5 2 6 3 4 5 3];
hold on; plot(x(I),y(I),'-'); hold off;
axis equal off;
subplot(1,2,2);
x = x(1:end-1); y = y(1:end-1);
plot(x,y,'o','MarkerSize',8);
t = delaunay(x,y);
hold on; triplot(t,x,y); hold off;
axis equal off;
if doprint
   print -deps2 tri-maxmin;
end

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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
