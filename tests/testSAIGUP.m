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

% testTopSurfaceGrid
%close all

gravity([0 0 9.8])
gravity on;
test_grid=false;
int_vert=false;
int_vert_poro=false;
semi_implicit=false;
injected_volume=30e6;h0=40;
rate=1e6;
%100*milli*darcy*10*(diff(fluid.rho))/fluid.mu(1)*year*30*0.2/0.3
if test_grid
   nx=40;ny=120;nz=1;
   nx=10;ny=10;nz=2;
   dims= [nx,ny,nz];
   %physdims=dims.*[3000,9000,80];
   physdims=[3000,9000,80];
   g = cartGrid(dims,physdims); %, [1 1 1]);
   %g.nodes.coords(:,[1,2])=bsxfun(@minus,g.nodes.coords(:,[1,2]),physdims(:,[1,2])/2);
   r2=(g.nodes.coords(:,1).^2+g.nodes.coords(:,2).^2);
   %g.nodes.coords(:,3)=g.nodes.coords(:,3)-10*physdims(3).*(1-r2/max(r2));
   g.nodes.coords(:,3)=2000+g.nodes.coords(:,3)-0.2*g.nodes.coords(:,1);
   %g = removeCells(g, [1 4]);
   %g.nodes.coords = twister(g.nodes.coords);
   g = computeGeometry(g);
   %rock.perm = 100*rand(g.cells.num,1)*milli*darcy;
   rock.perm = ones(g.cells.num,1)*100*milli*darcy;
   rock.poro = 0.3*ones(g.cells.num,1);
   rock.poro(1:nx*ny) = 0.01
   %figure
   %plotGrid(g); view(3)
   grdecl = [];
else
   makeSAIGUP
   
   %load nograv
   %load norne
   grdecl = [];
end

g_top = topSurfaceGrid(g, 'grdecl', grdecl);
rock2D = averageRock(rock, g_top);

figure(1),clf;
subplot(1, 3, 1)
plotGrid(g_top)
subplot(1, 3, 2)
cn_top = cellNodes(g_top);
%plotCellData(g_top, accumarray(cn_top(:,1), g_top.nodes.z(cn_top(:,3)))./accumarray(cn_top(:,1),1));
%plotCellData(g_top, g_top.cells.z)
colorbar
ca = caxis;

subplot(1, 3, 3)
plotCellData(g, g.cells.centroids(:,3))
view(3);
colorbar
%return
% 
% figure;
% plotCellData(g_top, rock2D.perm)
% 
H = max(g_top.cells.H);
fluid = initVEFluid(g_top, 'mu', [0.1 0.4]*centi*poise, ...
         'rho', [600 1000],'sr', 0.1,'sw',0.0);
%fluid = initVEFluid(g_top);
sol = initResSol(g_top, 0);
sol.h = zeros(g_top.cells.num, 1);
total_time= 30*year;
dt = total_time/50;
%find(min(g_top.cells.z));
%sol.h(g_top.cells.z==max(g_top.cells.z)) = min(g_top.cells.H)/10;
%ind_well=cart2active(g_top,sub2ind(g_top.cartDims,6,60));
well_pos=[0.412,4.46]*1e3;
%well_pos=[1.5,4.5]*1e3;
r2=sum(bsxfun(@minus,g_top.cells.centroids,well_pos).^2,2);
[mm,ind_well]=min(r2);
%well_pos=g_top.cells.centroids(ind_well,:);
sigma=10*injected_volume/(pi*h0);%assume average porosity 0.1
sol.h= h0*exp(-r2/sigma);
sol.h=sol.h*0.0;
sol.max_h=sol.h;

%sol.h(g_top.cells.num/2:(g_top.cells.num/2)+10) = min(g_top.cells.H)/10;
%sol.h(g_top.cells.num/2+5) = min(g_top.cells.H)/10;
%sol.h(end) = min(g_top.cells.H)/10;
%sol.h(1) = min(g_top.cells.H)/10;
t = 0;

bc = pside([], g_top, 'LEFT', 100*barsa())
bc = rmfield(bc,'sat');
bc.h=zeros(size(bc.face));
%bc =[];
W = addWell([], g_top, rock2D, ind_well,'Type', 'rate','Val',rate/year,'Radius',0.1)
%W = addWell(W, g_top, rock2D, g_top.cells.num,'Type', 'rate','Val',-rate/year,'Radius',0.1);bc=[];
%W=[];
for i=1:numel(W)
   %W(i) = rmfield(W(i),'compi');
   W(i).compi=nan;
   W(i).h = g_top.cells.H(W(i).cells);
end



%src = addSource(src, cells, values);
%set(gcf,'WindowStyle','docked')
figure(12);
%subplot(1, 2, 1)
%plotCellData(g_top, sol.h)
colorbar
%T_2d=computeTrans(g_top,rock2D)
S_2d = computeMimeticIPVE(g_top, rock2D,'Innerproduct','ip_simple');
preComp = initTransportVE(g_top, rock2D); 

%S_2d_org = computeMimeticIP(g_top, rock2D,'Innerproduct','ip_tpf');
%return
while t < total_time
   % plot(sol.h)
   % hold off;
   % ylim([0 1])
   %
   % drawnow
   %
   % sum(sol.h*dx)
   %
   % find(sol.h>0.01)
   disp(['Time  ',num2str(t/year)])
   
   pv_3D(g_top.columns.cells)=rock.poro(g_top.columns.cells)...
      .*rldecode(g_top.cells.volumes,diff(g_top.cells.columnPos));
   
   %free_volume= sum(sol.h.*rock2D.poro.*g_top.cells.volumes);
   if(int_vert_poro)
      free_volume = sum(integrateVertically(pv_3D, sol.h, g_top));
   else
      free_volume= sum(sol.h.*rock2D.poro.*g_top.cells.volumes);
   end
   disp(['Total volume free ', num2str(free_volume/1e6)]);
   if(isfield(sol,'max_h'))
      if(int_vert_poro)
         trapped_volume= sum(integrateVertically(pv_3D, sol.max_h-sol.h, g_top))*fluid.sr;
      else
         trapped_volume=sum((sol.max_h-sol.h).*rock2D.poro.*g_top.cells.volumes)*fluid.sr;
      end
      disp(['Total volume trapped ', num2str(trapped_volume/1e6)]);
      disp(['Total volume ', num2str((trapped_volume+free_volume)/1e6)]);
   end
   gravity on
   sol = solveIncompFlowVE(sol, g_top, S_2d, rock, fluid, 'bc', bc,'wells',W);
   %%gravity on
   sol = explicitTransportVE(sol, g_top, dt, rock, fluid,...
      'Verbose',true,'intVert',int_vert,'intVert_poro',int_vert_poro, 'semi_implicit', semi_implicit,...
      'bc', bc,'wells',W,'preComp', preComp);
   figure(12);
   subplot(1, 2, 1)
   plotCellData(g_top, sol.h)
   caxis([0 max(sol.h)]);
   %colorbar
   subplot(1, 2, 2)
   plotCellData(g_top, sol.pressure(1:g_top.cells.num));
   caxis([50,150]*barsa); %colorbar
   % figure(10)
   % plot(g_top.cells.z, 'r*')
   % hold on;
   % plot(g_top.cells.z-sol.h)
   % hold off;
   %
   % figure(11)
   % plot(sol.h)
   % hold off;
   %
   
   figure(11)
   plotCellData(g, height2Sat(sol, g_top, fluid))
   view([-6 50])
   
   
   drawnow
   t = t + dt;
   
   
end


% 
% 
% figure;
% subplot(1, 3, 1)
% plotGrid(g_top)
% subplot(1, 3, 2)
% cn_top = cellNodes(g_top);
% plotCellData(g_top, accumarray(cn_top(:,1), g_top.nodes.z(cn_top(:,3)))./accumarray(cn_top(:,1),1));
% colorbar
% ca = caxis;
% 
% subplot(1, 3, 3)
% plotCellData(g, g.cells.centroids(:,3))
% colorbar
% caxis(ca);
% 
% 
% figure;
% 
% subplot(1, 2, 1)
% plotCellData(g_top, convertTo(rock2D.perm, milli*darcy))
% colorbar
% 
% subplot(1, 2, 2)
% plotCellData(g, convertTo(rock.perm(:,1), milli*darcy))
% colorbar
% view(3)
% 
% 



%g_top.nodes.z

% g_top.columns.cells
% g_top.columns.kPos
% g_top.columns.dz

%g_top.cells.normals
