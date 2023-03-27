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
close all
cases_made={'org'}
mycase=cases_made{1}
gravity([0 0 9.8])
gravity on;
test_grid=false;
int_vert=false;%true is not properly implemented for the pressuresolver
int_vert_poro=false;
semi_implicit=false;
injected_volume=4e6;h0=10;
rate=1e6;
% atestgrid for testing
if(~readCache({mycase},'verbose', true)) %#ok<NOPTS>
   switch mycase
      case 'test_case'
         %nx=120;ny=200;nz=1;
         nx=5;ny=5;nz=1;
         dims= [nx,ny,nz];
         physdims=[6000,10000,25];
      g = cartGrid(dims,physdims); %, [1 1 1]);
      g.nodes.coords(:,2)=g.nodes.coords(:,2)-10000;
      %r2=(g.nodes.coords(:,1).^2+g.nodes.coords(:,2).^2);
      g.nodes.coords(:,3)=800+g.nodes.coords(:,3)-0.5*g.nodes.coords(:,1);
      g = computeGeometry(g);
      rock.perm = ones(g.cells.num,1)*100*milli*darcy;
      rock.poro = 0.3*ones(g.cells.num,1);
      grdecl = [];
      %    g_top = topSurfaceGrid(g, 'grdecl', grdecl);
      case 'org'
      %% check makeSome2DGrid for manipulations
      % now the two upper layers are removed
      [g_top,g,rock,grdecl]=makeSome2DGrids('sleipner');
      %% if one will experiment with flat surface
      %g_top.cells.z=700*ones(size(g_top.cells.z));
      %g_top.faces.z=700*ones(size(g_top.faces.z));
      %g_top.nodes.z=700*ones(size(g_top.faces.z));
      %rock2D = averageRock(rock, g_top);
      %g_top.nodes.z=700*ones(size(g_top.nodes.z));
      %% if one will use constant permeability
      %rock.perm = ones(g.cells.num,1)*mean(rock2D.perm);
      %rock.poro = ones(g.cells.num,1)*mean(rock2D.poro);
      %rock2D = averageRock(rock, g_top);
      clear g;
      %load norne
      otherwise
         error('not valid')
   end

   
   writeCache({mycase}, {'g_top','rock','grdecl'});
end
%%
% define the cartesian coordinates for nicer plotting
X=reshape(g_top.cells.centroids(:,1),g_top.cartDims(1),g_top.cartDims(2));
Y=reshape(g_top.cells.centroids(:,2),g_top.cartDims(1),g_top.cartDims(2));
Z=reshape(g_top.cells.z,g_top.cartDims(1),g_top.cartDims(2));
rock2D = averageRock(rock, g_top);
PERM=reshape(rock2D.perm,g_top.cartDims(1),g_top.cartDims(2));
%return
%%
% time data for injection rate at the moment not used
a=xlsread('data/sleipner/Injection rates Layer 9.xls');
time_values=a(12:end,3);
injectiondata=a(12:end,4:6);
rates = diff(injectiondata(2:end,end));
density_surface=a(6,7);
density_reserviour=a(6,4);
%return
% 
figure;
plotCellData(g_top, rock2D.perm)
%%
%define fluid 
fluid = initVEFluid(g_top, 'mu', [0.1 0.4]*centi*poise, 'rho', [600 1000],'sr',0.0);
sol = initResSol(g_top, 0);
sol.h = zeros(g_top.cells.num, 1);
%% define time intervals (should be change to real data)
total_time= 30*year;
dt = total_time/(30*4);
%% I have calculated the position of the wells so that it is easier to
%% expermint with smaller grids. I find the well_cell in the the other
%% makeSome2dGrids can be used to slice the model by differnet parameterst
%% to cutGRDECL
%find(min(g_top.cells.z));
%sol.h(g_top.cells.z==max(g_top.cells.z)) = min(g_top.cells.H)/10;
%I=61; J=126, K=7
%ind_well_org=cart2active(g_top,sub2ind(g_top.cartDims,61,126));
%well_pos=g_top.cells.centroids(ind_well,1:2);
well_pos= [3.031296727824156   -3.722750260389535]*1e3;
r2=sum(bsxfun(@minus,g_top.cells.centroids,well_pos).^2,2);
[mm,ind_well]=min(r2);
%well_pos=g_top.cells.centroids(ind_well,:);
%% if no injection rate is present I put some CO2 in the system
sigma=10*injected_volume/(pi*h0);%assume average porosity 0.1
sol.h= h0*exp(-r2/sigma);
if(rate>0)
   sol.h=sol.h*0.0;
end
%% h is always smaller than H
sol.h=min([sol.h,g_top.cells.H],[],2);
sol.max_h=sol.h;
t = 0;
%%set hydrostatic boundary conditions on all sided
% is this good ?
bc_faces=find(sum(g_top.faces.neighbors>0,2)==1);
bc = addBC([], bc_faces, 'pressure', g_top.faces.z(bc_faces)*1000*norm(gravity));
bc = rmfield(bc,'sat');
bc.h=zeros(size(bc.face));
W = addWell([], g_top, rock2D, ind_well,'Type', 'rate','Val',rate/year,'Radius',0.1);
%% use MRST to add wells we need to change defintion of well
for i=1:numel(W)
   W(i).compi=nan;
   W(i).h = g_top.cells.H(W(i).cells);
end
figure(10);
colorbar
%one can choose between different innerporduct for the mimetic calculation
S_2d = computeMimeticIPVE(g_top, rock2D,'Innerproduct','ip_simple');
%S_2d = computeMimeticIP(g_top, rock2D,'Innerproduct','ip_tpf');
%% time loop start with simple time control
while t < total_time
   disp(['Time  ',num2str(t/year)])
   cum_z  = cumulativeHeight(g_top);
  
   if(int_vert_poro)
      pv_3D(g_top.columns.cells)=rock.poro(g_top.columns.cells)...
         .*rldecode(g_top.cells.volumes,diff(g_top.cells.columnPos));
      free_volume = sum(integrateVertically(pv_3D, sol.h, g_top, cum_z));
   else
      free_volume= sum(sol.h.*rock2D.poro.*g_top.cells.volumes);
   end
   disp(['Total volume free ', num2str(free_volume/1e6)]);
   if(isfield(sol,'max_h'))
      if(int_vert_poro)
         trapped_volume= sum(integrateVertically(pv_3D, sol.max_h-sol.h, g_top, cum_z))*fluid.sr;
      else
         trapped_volume=sum((sol.max_h-sol.h).*rock2D.poro.*g_top.cells.volumes)*fluid.sr;
      end
      disp(['Total volume trapped ', num2str(trapped_volume/1e6)]);
      disp(['Total volume ', num2str((trapped_volume+free_volume)/1e6)]);
   end
   sol = solveIncompFlowVE(sol, g_top, S_2d, rock, fluid, 'bc', bc,'wells',W);

   %% some plotting for loocking at initial pressure
   if false %(true)
      %%
      cellflux = faceFlux2cellFlux(g_top,sol.flux);
      cellvel  = cellFlux2cellVelocity(g_top,cellflux);
      vel=bsxfun(@times,g_top.faces.normals,sol.flux);
      vel=bsxfun(@rdivide,vel,sqrt(sum(vel.^2,2)));
      % figure(22)
      % pos=g_top.faces.centroids;
      % clf,plotGrid(g_top),hold on,quiver(pos(:,1),pos(:,2),vel(:,1),vel(:,2),0.5)
      % some velocity plotting
      figure(33)
      plotVelocityVE(g_top,sol)
      % 
      figure(44)
      vx=reshape(cellvel(:,1),g_top.cartDims);vy=reshape(cellvel(:,2),g_top.cartDims);
      v_norm=sqrt(vx.^2+vy.^2);
      %startcells=find(sqrt(r2)<100);
      [I,J]=ndgrid(1:floor(g_top.cartDims(1)/30):g_top.cartDims(1),...
         1:floor(g_top.cartDims(2)/6):g_top.cartDims(2));
      startcells = sub2ind(g_top.cartDims,I(:),J(:));
      pos=g_top.cells.centroids(startcells,:);
      %r=100;vv=0:01:pi;clf,streamline(X',Y',(vx./v_norm)',(vy./v_norm)',pos(:,1),pos(:,2))
      cnz=g_top.cells.normals(:,3)./sqrt(sum(g_top.cells.normals.^2,2));
      %plot situation as it you see a landscape with lakes.
      figure(99),
      surf(X,Y,Z+reshape(sol.h.*cnz,g_top.cartDims),reshape(sol.h,g_top.cartDims))
      set(gca,'zdir','reverse');view([1 1 0.0]);shading interp
      c=stream2(X',Y',(vx./v_norm)',(vy./v_norm)',pos(:,1),pos(:,2));
      b=stream2(X',Y',-(vx./v_norm)',-(vy./v_norm)',pos(:,1),pos(:,2));
      a ={c{:},b{:}};
      xy2d=vertcat(a{:});     
      if(false)
         zz=griddata(X',Y',Z',xy2d(:,1),xy2d(:,2));
         %%
         pos_start=1;
         for kk=1:numel(a)
            pos_end=pos_start+size(a{kk},1)-1;
            aa{kk}=[a{kk},zz(pos_start:pos_end)];
            pos_start=pos_end+1;
         end
         %%
         figure(44),clf
         clf,
         mesh(X,Y,Z,reshape(sol.h,g_top.cartDims),'EdgeAlpha',0.25,'FaceAlpha',0.375)
         %mesh(X,Y,Z,sqrt(vx.^2+vy.^2))
         view([1 1 0.0]);set(gca,'zdir','reverse')%shading interp         
         hold on,streamline(aa),grid off
      end
      %alpha(0.1)
      %%
      
      return
   end
   %%Plot landscape with water C02 lakes
   cnz=g_top.cells.normals(:,3)./sqrt(sum(g_top.cells.normals.^2,2));
   figure(99),
   if(false)
      clf,surf(X,Y,Z+reshape(sol.h.*cnz,g_top.cartDims),reshape(sol.h,g_top.cartDims))
      set(gca,'zdir','reverse');view([1 1 0.0]);shading interp
   else
      clf,pcolor(X,Y,reshape(sol.h.*cnz,g_top.cartDims))
   end
   %%
   %gravity on
   sol = explicitTransportVE(sol, g_top, dt, rock, fluid,...
      'Verbose',true,'intVert',int_vert,'intVert_poro',int_vert_poro, 'semi_implicit', semi_implicit,...
      'bc', bc,'wells',W);
   figure(10);
   subplot(3, 2, 1)
   mesh(X,Y,-Z,reshape(sol.h,g_top.cartDims));colorbar
   subplot(3, 2, 2)
   mesh(X,Y,reshape(sol.h,g_top.cartDims));colorbar
   %plotCellData(g_top, sol.h)
   caxis([0 max(g_top.cells.H)]);
   %colorbar
   subplot(3, 1, 3)
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
   
   %%
   %figure(11),cla
   %plotGrid(g,'FaceColor','none','Edgealpha',0.05);
   %clf,cla,plotGrid(g,find(height2Sat(sol, g_top, fluid)>0.2),'FaceColor','none','Edgealpha',0.05);
   %plotCellData(g, height2Sat(sol, g_top, fluid),find(height2Sat(sol, g_top, fluid)>0.1))
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
