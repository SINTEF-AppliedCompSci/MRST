doprint = true;

%% Single cell
c1 = [0.9 0.9 0.9];
c2 = [0.8 0.8 0.8];
d = 0.2;
clf;
X = [0 1 1 0; 0 1 1 0; 0 0 0 0; 1 1 1 1; 0 1 1 0; 0 1 1 0]';
Y = [0 0 1 1; 0 0 1 1; 0 0 1 1; 0 0 1 1; 0 0 0 0; 1 1 1 1]';
Z = [0 0 0 0; 1 1 1 1; 0 1 1 0; 0 1 1 0; 0 0 1 1; 0 0 1 1]';
patch(X,Y,Z,c1);
view(-70,30); 
axis([-0.05 1.05 -0.05 1.05 -0.05 1.05]+d*repmat([-1 1],1, 3));
axis off;
if doprint
   print -deps2 cell-def;
end

%% Cell with faces
Xp = X; Xp(:,3) = Xp(:,3)-d; Xp(:,4)=Xp(:,4)+d;
Yp = Y; Yp(:,5) = Yp(:,5)-d; Yp(:,6)=Yp(:,6)+d;
Zp = Z; Zp(:,1) = Zp(:,1)-d; Zp(:,2)=Zp(:,2)+d;
patch(Xp,Yp,Zp,c2);
hold on;
plot3([X(:) Xp(:) NaN*ones(24,1)]', [Y(:) Yp(:) NaN*ones(24,1)]', ...
   [Z(:) Zp(:) NaN*ones(24,1)]', ':k');
hold off;
if doprint
   print -deps2 faces-def;
end

%% Cell with edges and vertices
clf;
plot3(X,Y,Z,'ok','LineWidth',2,'MarkerSize',8);
hold on;
plot3(X([1:end 1],:),Y([1:end 1],:),Z([1:end 1],:),'-k');
hold off;
view(-70,30); 
axis([-0.05 1.05 -0.05 1.05 -0.05 1.05]+d*repmat([-1 1],1, 3));
axis off;
if doprint
   print -deps2 vertices-def;
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
