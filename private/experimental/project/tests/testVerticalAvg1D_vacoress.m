
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

gravity on

   
%% Make grid   
nn=10; %100;
L=5000;
%volume=10e6;
volume = 0.01e6; h0=1; 
%h0=15; 
H=15;
n_well=nn+1;

dx=L/(2*nn+1);
dy=1000;




theta=10*pi/180; %1*pi/180;
x_grid=[[n_well-1:-1:1]*(-dx),0,[1:n_well-1]*dx];
N=length(x_grid);
z_grid=sin(theta).*x_grid;
z_grid(:) = z_grid(end:-1:1);   
z = [(z_grid-H), (z_grid-H), z_grid, z_grid]';        
g = cartGrid([N-1, 1, 1], [L-L/N dy H]); %, [1 1 1]);
g.nodes.coords(:,3) =  z;
%  g = removeCells(g, [1 4]);


g = computeGeometry(g);
grdecl = [];
g_top = topSurfaceGrid(g, 'grdecl', grdecl);
%% rock and fluid data
% rock.perm = 100*rand(g.cells.num,1)*milli*darcy;
% rock.perm = (1:g.cells.num)'*milli*darcy;
K=(100e-3)*darcy();
phi= 0.1;

%K = 1; phi = 1;

%n_well=100;
rock.perm = K*ones(g.cells.num,1);
rock.poro = phi*ones(g.cells.num,1);

fluid = initVEFluid(g_top, 'mu', [0.1 0.4]*centi*poise, 'rho', [600 1000], 'sr', 0, 'sw', 0);
sol = initResSol(g_top, 0);


sigma=dy*h0*sqrt(pi)/(volume/phi);
h=zeros(N,1)+min(h0*exp(-(sigma*(x_grid+L/3))'.^2),H);
sol.h = h(1:end-1);



figure(1);
total_time=5*year;
dt = total_time/1;
%set(gcf,'WindowStyle','docked')
t=0;
sol.flux(1:g_top.cells.num+1)=1e-30;
bc=[];
bc=addBC(bc,1,'pressure',1,'sat',0);
bc=addBC(bc,g_top.cells.num+1,'pressure',1,'sat',0);
for i=1:numel(bc)
   bc(i).h=0;
end
vacoress=true;
numthreads=4;
if(vacoress)
   addpath('/home/hnil/heim/SVN/VACORESS/vacoress_new/build/mex');
   %addpath('/home/hnil/heim/SVN/vacoress_mine/vacoress/build/mex')
   VETransportCPU();
   %fastTransportVE()
   sol.max_h=sol.h;
end
sol_old=sol;

while t < total_time
   tic;
   if(~vacoress)
      sol = explicitTransportVE(sol, g_top, dt, rock, fluid,'time_stepping','coats','Verbose',false,'bc',bc);
   else
      %{
    [sol.h,sol.max_h] = fastTransportVE(sol, g_top, dt, rock, fluid, 'computeDt',true,...
        'bc', [], 'gravity', norm(gravity), 'flag', numthreads);
      %}
    [sol.h,sol.max_h] = VETransportCPU(sol, g_top, dt, rock, fluid, 'computedt',...
        false,'verbose',true,'gravity', norm(gravity), 'flag', numthreads); 
   end
   toc;
   plot(sol.h)
   hold on
   plot(sol.max_h, 'r');
   %plot(sol.max_h-sol.h, 'g');
   hold off;

   sum(sol.h.*rock.poro);
   %sum(sol.h.*rock.poro.*(1-fluid.sw)+rock.poro.*(sol.max_h-sol.h)*fluid.sr)
   %ylim([0 1])

   drawnow

   % figure(11)
   % plot(sol.h)
   % hold off;
   %
   % drawnow

   t = t + dt;

end


