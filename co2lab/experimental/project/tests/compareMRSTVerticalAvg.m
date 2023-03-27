% testTopSurfaceGrid
clear, close all
gravity on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% Quantities to vary:                        %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 30;    % number of cells in x-direction
L = 1000;  % physical dimension in x-direction
dy = 500; % physical dimension in y-direction
volume = 0.05e6; %volume of CO2 plume
h0=8; % max height of initial CO2 plume

total_time  = 40*year; 
% need small timesteps to avoid instabilities in 3D solution because the
% plume is thin
dt = year/10;

% flux-function, 1 = linear flux function
n_flux = 1;

% set resolution in z-direction
dim3 = 10;

% refine near top for the 3D grid in the z direction
refine = true; 

% height of reservoir
H = 15;
% dip angle of reservoir
theta = 5*pi/180; 


%% Make grid   
% We first make a regular grid
dx     = L/(n+1); 
x_grid = [(n/2:-1:1)*(-dx),0,(1:n/2)*dx];

% define variation in z-coordinates according to theta
z_grid = sin(theta).*x_grid;
z_grid(:) = z_grid(end:-1:1);   


g = cartGrid([n, 1, dim3], [L dy H]);
g = computeGeometry(g);

% number of nodes in one layer
numN = (g.cartDims(1)+1)*(g.cartDims(2)+1);
%% Refine 3D grid near the top. 
% The 3D grid must be quite fine to capture the CO2 plume. For efficiency
% we refine more in the upper part of the grid we have CO2.
if refine
   for k = 2:g.cartDims(3)
      ix = ((k-1)*numN+1):(k*numN);
      g.nodes.coords(ix,3) = g.nodes.coords(ix,3) - ...
                         (H-g.nodes.coords(ix,3))*(k/(g.cartDims(3)-1)).^2;
   end
   figure(1);
   g = computeGeometry(g);
  
end
g_plot = g; g_plot = computeGeometry(g_plot);

g.nodes.coords(:,3) = g.nodes.coords(:,3) + repmat(z_grid', 2*(dim3+1), 1);

g = computeGeometry(g);

maxx=max(g.nodes.coords(:,1));

g_top = topSurfaceGrid(g);

%% Set rock and fluid data
%K = rand(dim3,1)*10+100;
%K = bsxfun(@times, ones(g.cartDims(1), g.cartDims(3)), K')*milli*darcy; logNormLayers(g.cartDims, dim3/2)*milli*darcy;
K = (100e-3)*darcy(); 
phi = 0.1;

%rock.perm  = reshape(K, [],1); %*ones(g.cells.num,1);
rock.perm  = K*ones(g.cells.num,1);
rock.poro  = phi*ones(g.cells.num,1);
fluid      = initVEFluid(g_top, 'mu', [0.1 0.4]*centi*poise, ...
                         'rho', [600 1000],'sr',0.0,'sw',0.0);
fluid_mrst = initSimpleFluid('mu', [ 0.1,  0.4]*centi*poise, ...
                             'rho', [600 1000], 'n', [n_flux n_flux]);

drho = 400;

%% Show 3D grid
clf
plotGrid(g); view(3), axis tight
xlabel('x'), ylabel('y'), zlabel('z')
view([-14 16])

%% Show 2D/topsurface grid
clf
plotGrid(g_top, 'faceColor', 'b'); 
view([-14 16]), axis tight

xlabel('x'), ylabel('y'), zlabel('z')
%% Display permeability field on 3D grid with theta = 0,
% NB: observe scaling of axis
clf
plotGrid(g_plot); view([0 0])
%plotCellData(g_plot, rock.perm); view([0 0])
xlabel('x'), ylabel('y'), zlabel('z')


%% Initialize solution structure
sol = initResSol(g_top, 0);
sol_mrst = initResSol(g, 0);

% initial CO2 plume
sigma=dy*h0*sqrt(pi)/(volume/phi);
h = zeros(n+1,1)+min(h0*exp(-(sigma*(x_grid+L/3))'.^2),H);

sol.h = h(1:end-1);
sol.s = height2Sat(sol, g_top, fluid);

sol_mrst.s = sol.s;
sol_mrst.h = sol.h;

% compute inner product
S = computeMimeticIP(g, rock);


%% Run transport-simulation:
figure(1);
error_s = 0;
error_h = 0;
xyz(1, :)      = findMassCenter(g, sol); 
xyz_mrst(1, :) = findMassCenter(g, sol_mrst); 
t = 0;
i = 2;
while t < total_time
 
subplot(3, 1, 1)
cells = 1:g_top.cells.num;
plot(cells, sol.h, cells , sol_mrst.h)
legend('2D', '3D');
ylim([0 h0*1.3]);


subplot(3, 1, 2)
plotCellData(g_plot, height2Sat(sol, g_top, fluid));
view([0 0]); caxis([0 1]);
axis tight off;
%xlim([maxx-3*dx maxx])
title(['vertAvg saturation ', 't = ', num2str(t/year),  'years']); 

subplot(3, 1, 3)
%sat_mrst=sol_mrst.h(:,1);
plotCellData(g_plot, sol_mrst.s(:,1));
axis tight off;
%xlim([maxx-3*dx maxx])
view([0 0]); caxis([0 1]); 
title(['mrst saturation']);
%zlim([H-h0 H])
drawnow  

sol_mrst = solveIncompFlow(sol_mrst, g, S, fluid_mrst);
sol_mrst = explicitTransport(sol_mrst, g, dt, rock, fluid_mrst, ...
                             'max_dt', dt, 'verbose', false);

sol = explicitTransportVE(sol, g_top, dt, rock, fluid, ...
                               'computeDt', true, 'intVert_poro', false);

% convert
sol.s = height2Sat(sol, g_top, fluid);
sol_mrst.h = sat2height(sol_mrst.s, g_top, rock);

xyz_mrst(i, :) = findMassCenter(g, sol_mrst); 
xyz(i, :)      = findMassCenter(g, sol); 


error_s = [error_s; norm(sol.s-sol_mrst.s(:,1))./norm(sol_mrst.s(:,1))];
error_h = [error_h; norm(sol.h-sol_mrst.h)./norm(sol.h)];
t = t + dt; 
i = i+1;


%  if(isStable(sol, g_top)), disp(['vertical stable, t = ', num2str(t/year)]), end
%  if(isStable(sol_mrst, g_top)),disp(['mrst stable, t = ', num2str(t/year)]); break; end
%  
 

end

time = 0:dt/year:total_time/year;

figure(1)
subplot(2,1,1)
plot(time, error_s,time, error_h); legend('e_s', 'e_h');
title(['dx = ', num2str(round(dx)),', Dt = ', num2str(dt/year), '. Dim3 = ', num2str(dim3)]); 
subplot(2, 1, 2)
plot(time(2:end), diff(xyz(:,1)), 'b', time(2:end), -diff(xyz(:,3)), 'g',  ...
     time(2:end), diff(xyz_mrst(:,1)), 'r', time(2:end), -diff(xyz_mrst(:,3)), 'k');
legend('vx', 'vz', 'vx-mrst', 'vz-mrst');   
title('Mass center velocity'); 


% %% compute timestep
% sinth= (z_grid(2:end)-z_grid(1:end-1))'/dx; sinth= [sinth(1);sinth;sinth(end)];
% const_perp = norm(gravity())*K*drho*sqrt(1-sinth.^2);
% const_par = norm(gravity())*K*drho*sinth;
% t_inject=20*year; t = 0; hh=0:H/100:H;
% v=[-well_rate*ones(n_well,1);well_rate*ones(n_well,1)]/(2*dy*H); c2 = 0;
% lam1 = @(h,h_res) h/fluid.mu(1); lam2 = @(h,h_res) ((H-h)+h_res*c2)/fluid.mu(2);
% FF = @(h,h_res) lam1(h,h_res).*lam2(h,h_res)./(lam1(h,h_res)+lam2(h,h_res));
% dd_max=max(abs( (FF(hh(1:end-1),0)-FF(hh(2:end),0))./(hh(1:end-1)-hh(2:end))));
% dt_perp=min(0.1./(2*pi*abs(const_perp)*dd_max*max(H)./(dx^2)));
% dt_par = min(0.5./(abs(const_par)*dd_max/(dx*phi)));
% dt_adv=min(0.1./(well_rate/(dx*dy*H)));
% %NB: 
% %dt = abs(min([dt_perp,dt_par,dt_adv]))
% 



