X = [-6 -5.6 -5.2 -5 -4.5  -4   -3.7 -3.2 -2.7 -2   -1.5  1.5    2  2.5    3  3.5   4];
Y = [1 0.52 0.65 0.21 0.57 0.11 0.64 0.28 0.74 0.25 0.35 0.35 0.31 0.39 0.31 0.44 0.32]; 
x=-6:0.05:4;
plot(10.^x, interp1(X([1:11 end]),Y(1:12),x,'pchip'), '--r', ...
   10.^x, interp1(X,Y,x,'pchip'), '-b','LineWidth', 2);
set(gca,'XScale','log', 'YLim', [0 1]);
h=legend('homogeneous medium','inhomogeneous medium',4);
set(h,'FontSize',12);
hold on, 
plot([0.03 0.03],[0 1],'--k','LineWidth',1);
plot([40 40],[0 1],'--k','LineWidth',1);
text(0.3,0.8,'REV','FontSize',16);
text(1e-5,0.85,'Microscopic effects','FontSize',12);
hold off

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
