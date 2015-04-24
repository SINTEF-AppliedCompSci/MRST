%{
Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.

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
%%
[g_top,g,rock,grdecl]=makeSome2DGrids('sleipner');
%%
X=reshape(g_top.cells.centroids(:,1),g_top.cartDims(1),g_top.cartDims(2));
Y=reshape(g_top.cells.centroids(:,2),g_top.cartDims(1),g_top.cartDims(2));
Z=reshape(g_top.cells.z,g_top.cartDims(1),g_top.cartDims(2));
rock2D = averageRock(rock, g_top);
PERM=reshape(rock2D.perm,g_top.cartDims(1),g_top.cartDims(2));
%return
a=xlsread('data/sleipner/Injection rates Layer 9.xls')
times=a(11:end,3);
injectiondata=a(11:end,4:6);
density_surface=a(6,7);
density_reserviour=a(6,4);
plot(3)
subplot(3,1,1)
plot(times,injectiondata(:,3))
subplot(3,1,2)
plot(diff(diff(injectiondata(:,3))))
subplot(3,1,3)%of some strange reson the 3 derivative of commulative colume is constant
plot(diff(diff(diff(injectiondata(:,3)))),'*')
%%
%plot in figure 1 and 2
plotRelPerm(g_top,rock)
figure(4)%the permeability is almost only height correlated and the z permeability is 10 times permX
plot(g.cells.centroids(:,3),rock.perm(:,3),'.')
stats = @(v) [min(v); max(v); mean(v); std(v)]
orgperm=max(rock.perm(:,1))...
   -(max(rock.perm(:,1))-min(rock.perm(:,1)))/(max(g.cells.centroids(:,3))-min(g.cells.centroids(:,3)))...
   *(g.cells.centroids(:,3)-min(g.cells.centroids(:,3)));
disp('Statistics of perm')
stats(convertTo(rock.perm, milli*darcy))
stats(convertTo(rock.perm(:,1)-orgperm, milli*darcy))
figure(5)
mesh(X,Y,Z,PERM)
