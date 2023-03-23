function drawAxisCross(s)
%Undocumented function

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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

if nargin<1, s=.2; end
a = axis;
[xm,ym,zm] = deal(a(1),a(3),a(5));
[xM,yM,zM] = deal(a(2),a(4),a(6));
e = ones(3,1);
v = [s*(xM-xm), s*(yM-ym), s*(zM-zm)];
V = diag(v);
hold on;
quiver3(xm-v(1)*e,ym-v(2)*e,zm-0*e,V(:,1),V(:,2),V(:,3),...
   'LineWidth',2,'Color','k');
text(xm-v(1)*e(1)+V(1,1), ym-v(2)*e(1)+V(1,2), zm+V(1,3), 'x','FontSize',12);
text(xm-v(1)*e(1)+V(2,1), ym-v(2)*e(1)+V(2,2), zm+V(2,3), 'y','FontSize',12);
text(xm-v(1)*e(1)+V(3,1), ym-v(2)*e(1)+V(3,2), zm+V(3,3), 'z','FontSize',12);
hold off
end
