
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

%n_well=100;
well_rate=0.0*1e6/year;
nn=100; %100;
L=5000;
%volume=10e6;
volume = 0.01e6; h0=1; 
%h0=15; 
H=15;
n_well=nn+1;

dx=L/(2*nn+1);
dy=1000;




theta=1*pi/180; %1*pi/180;

x_grid=[[n_well-1:-1:1]*(-dx),0,[1:n_well-1]*dx];
N=length(x_grid);
z_grid=sin(theta).*x_grid;
%z_grid2=z_grid-sin(theta)*L*(x_grid-L/2).*(x_grid+L/2)/(L^2/4);
z_grid(:) = z_grid(end:-1:1);   
z = [(z_grid-H), (z_grid-H), z_grid, z_grid]';     
   
g = cartGrid([N-1, 1, 1], [L-L/N dy H]); %, [1 1 1]);

g.nodes.coords(:,3) =  z;
%  g = removeCells(g, [1 4]);


g = computeGeometry(g);
g_top = topSurfaceGrid(g);
%% rock and fluid data
% rock.perm = 100*rand(g.cells.num,1)*milli*darcy;
% rock.perm = (1:g.cells.num)'*milli*darcy;
K=(100e-3)*darcy();
phi= 0.1;

%K = 1; phi = 1;

%n_well=100;
rock.perm = K*rand(g.cells.num,1);
rock.poro = phi*rand(g.cells.num,1);

fluid = initVEFluid(g_top, 'mu', [0.1 0.4]*centi*poise, 'rho', [600 1000], 'sr', 0, 'sw', 0);
drho=400;
%lam1 = @(h) h/mu1; lam2 = @(h) (H-h)/mu2;
%FF = @(h) (1/(mu1*mu2))*(h.*(H-h))./(h/mu1 +(H-h)/mu2);
%FF = lam1;


sol = initResSol(g_top, 0);

%source(n_well)=0*2e6/year();
sigma=dy*h0*sqrt(pi)/(volume/phi);
%h=min(h0*exp(-(sigma*x_grid)'.^2),H);
h=zeros(N,1)+min(h0*exp(-(sigma*(x_grid+L/3))'.^2),H);
sol.h = h(1:end-1);


%zeros(N,1);h(n_well)=h0;
%const_perp = norm(gravity())*K*drho*cos(theta);
%const_par = norm(gravity())*K*drho*sin(theta);
sinth= (z_grid(2:end)-z_grid(1:end-1))'/dx;
sinth= [sinth(1);sinth;sinth(end)];
const_perp = norm(gravity())*K*drho*sqrt(1-sinth.^2);
const_par = norm(gravity())*K*drho*sinth;
tic;
%hh=0:H/100:H;

%dt=1000*year();% 1000*year()/10;
total_time=500*year;

dt = total_time/10;
t_inject=20*year;
t = 0;
v=[-well_rate*ones(n_well,1);well_rate*ones(n_well,1)]/(2*dy*H);
c2 = 0;
lam1 = @(h,h_res) h/fluid.mu(1);
lam2 = @(h,h_res) ((H-h)+h_res*c2)/fluid.mu(2);
FF = @(h,h_res) lam1(h,h_res).*lam2(h,h_res)./(lam1(h,h_res)+lam2(h,h_res));
hh=0:H/100:H;
dd_max=max(abs( (FF(hh(1:end-1),0)-FF(hh(2:end),0))./(hh(1:end-1)-hh(2:end))));
dt_perp=min(0.1./(2*pi*abs(const_perp)*dd_max*max(H)./(dx^2)));
dt_par = min(0.5./(abs(const_par)*dd_max/(dx*phi)));
dt_adv=min(0.1./(well_rate/(dx*dy*H)));
%abs(min([dt_perp,dt_par,dt_adv]))

%NB: 
dt=abs(min([dt_perp,dt_par,dt_adv]));



% 
% figure;
% subplot(1, 2, 1)
% plotGrid(g)
% view(3)
% subplot(1, 2, 2)
% plotGrid(g_top); 


figure;
%set(gcf,'WindowStyle','docked')
while t < total_time
   
   




% sum(sol.h*dx)
% 
% find(sol.h>0.01)


sol = explicitTransportVE(sol, g_top, dt, rock, fluid, 'computeDt', true,'time_stepping', 'dynamic');


% subplot(1, 2, 1)
% plotCellData(g_top, sol.h)
% colorbar
% subplot(1, 2, 2)

% figure(10)
% plot(g_top.cells.z, 'r*')
% hold on;
% plot(g_top.cells.z-sol.h)
% hold off;
% 
plot(sol.h)
hold on
plot(sol.max_h, 'r');
plot(sol.max_h-sol.h, 'g');
hold off;
legend({'h', 'hmax', 'hmax-h'})
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

