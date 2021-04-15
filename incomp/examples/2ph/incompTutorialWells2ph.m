%% Basic Transport-Solver Tutorial
% Consider a two-phase oil-water problem. Solve the two-phase pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\lambda K\nabla p,$$
%
% where v is the Darcy velocity (total velocity) and lambda is the
% mobility, which depends on the water saturation S.
%
% The saturation equation (conservation of the water phase) is given as:
%
% $$ \phi \frac{\partial S}{\partial t} +
%     \nabla \cdot (f_w(S) v) = q_w$$
%
% where phi is the rock porosity, f is the Buckley-Leverett fractional
% flow function, and q_w is the water source term.
%
% This is a continuation of the <matlab:edit('incompIntro.m') incompIntro>
% flow-solver tutorial</a>, in which we solved the corresponding
% single-phase problem using the two-point pressure solver. Here, we
% demonstrate how this flow solver can be extended by an explicit or an
% implicit two-phase transport solver. The grid is Cartesian with
% isotropic, homogeneous permeability. See the
% <matlab:edit('incompIntro.m') basic flow-solver tutorial</a> for more
% details on the grid structure, the structure used to hold the solutions,
% etc.

mrstModule add incomp

%% Define geometry and rock parameters
% Construct a Cartesian grid of size 20-by-20-by-5 cells, where each cell
% has dimension 1-by-1-by-1. Set the permeability K to be homogeneous,
% isotropic and equal 100 mD and the porosity to be equal to 0.3. Compute
% one-sided (half) transmissibilities from input grid and rock properties.
nx = 20; ny = 20; nz = 5;
G = cartGrid([nx, ny, nz]);
G = computeGeometry(G);

rock.perm  = 100*milli*darcy*ones(G.cells.num, 1);
rock.poro  = .3*ones(G.cells.num, 1);
rock.ntg   = ones(G.cells.num, 1);

hT  = computeTrans(G, rock, 'Verbose', true);

%% Define the two-phase fluid model
% The <matlab:help('initSimpleFluid') two-phase fluid model> has default values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 10] cP.
fluid = initSimpleFluid('mu' , [   1,  10] .* centi*poise     , ...
                        'rho', [1000, 700] .* kilogram/meter^3, ...
                        'n'  , [   2,   2]);

%%
% The fluid model represented by the <matlab:help('fluid') fluid structure>
% is the two-phase incompressible counterpart to the fluid models used in
% the <matlab:mrstExamples('ad-blackoil') black-oil framework>
%
s=linspace(0,1,20)'; kr=fluid.relperm(s);
plot(s, kr(:,1), 'b', s, kr(:,2), 'r', 'LineWidth',2);
title('Relative permeability curves')
legend('Water','Oil','Location','Best')

%% Introduce wells
% We will include two wells, one rate-controlled vertical well and one
% horizontal well controlled by bottom-hole pressure. Wells are described
% using a Peacemann model, giving an extra set of equations that need to be
% assembled, see the <matlab:edit('incompTutorialWells.m') tutorial on well
% models> for more details.
W = addWell([], G, rock, 1 : nx*ny : nx*ny*nz,          ...
            'InnerProduct', 'ip_tpf', ...
            'Type', 'rate', 'Val', 1.0/day(), ...
            'Radius', 0.1, 'Comp_i', [1, 0]);
W = addWell(W, G, rock, nx : ny : nx*ny, ...
            'InnerProduct', 'ip_tpf', ...
            'Type', 'bhp' , 'Val', 1.0e5, ...
            'Radius', 0.1, 'Dir', 'y', 'Comp_i', [0, 1]);

% To check if the wells are placed as we wanted them, we plot them
clf
plotGrid(G, 'FaceColor', 'none'); view(3);
[ht, htxt, hs] = plotWell(G, W, 'radius', 0.1, 'height', 2);
set(htxt, 'FontSize', 16);

%%
% Once the wells are added, we can generate the components of the linear
% system corresponding to the two wells and initialize the solution
% structure (with correct bhp)
%
rSol = initState(G, W, 0, [0, 1]);

%% Solve initial pressure in reservoir
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir and the wells.
gravity off
rSol = incompTPFA(rSol, G, hT, fluid, 'wells', W);

% Report initial state of reservoir
subplot(2,1,1), cla
   plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa));
   title('Initial pressure'), view(3)

subplot(2,1,2), cla
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   plotCellData(G, accumarray(cellNo, ...
      abs(convertTo(faceFlux2cellFlux(G, rSol.flux), meter^3/day))));
   title('Initial flux intensity'), view(3)

%% Transport loop
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation. This procedure is repeated for a given number of time steps
% (here we use 15 equally spaced time steps). The error introduced by this
% splitting of flow and transport can be reduced by iterating each time
% step until e.g., the residual is below a certain user-prescribed
% threshold (this is not done herein).
T      = 300*day();
dT     = T/15;
dTplot = 100*day();  % plot only every 100th day
N      = fix(T/dTplot);
pv     = poreVolume(G,rock);

%%
% The transport equation will be solved by the single-point upstream method
% with either explicit or implicit time discretizations. Both schemes may
% use internal time steps to obtain a stable discretization. To represent
% the two solutions, we create new solution objects to be used by the
% solver with implicit transport step.
rISol = rSol;

%% Start the main loop
t  = 0; plotNo = 1; hi = 'Implicit: '; he = 'Explicit: ';
e = []; pi = []; pe = [];
while t < T,
   rSol  = explicitTransport(rSol , G, dT, rock, fluid, 'wells', W);
   rISol = implicitTransport(rISol, G, dT, rock, fluid, 'wells', W);

   % Check for inconsistent saturations
   s = [rSol.s(:,1); rISol.s(:,1)];
   assert(max(s) < 1+eps && min(s) > -eps);

   % Update solution of pressure equation.
   rSol  = incompTPFA(rSol , G, hT, fluid, 'wells', W);
   rISol = incompTPFA(rISol, G, hT, fluid, 'wells', W);

   % Measure water saturation in production cells in saturation
   e = [e; sum(abs(rSol.s(:,1) - rISol.s(:,1)).*pv)/sum(pv)]; %#ok
   pe = [pe; rSol.s(W(2).cells,1)' ];                 %#ok
   pi = [pi; rISol.s(W(2).cells,1)'];                 %#ok

   % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   if ( t < plotNo*dTplot && t <T), continue, end

   % Plot saturation
   heading = [num2str(convertTo(t,day)),  ' days'];
   r = 0.01;
   subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.48]), cla
   plotCellData(G, rSol.s(:,1));
   view(60,50), axis equal off, title([he heading])

   subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.48]), cla
   plotCellData(G, rISol.s(:,1));
   view(60,50), axis equal off, title([hi heading])
   drawnow

   plotNo = plotNo+1;
end

%%
% As we clearly can see from the plots in the figure, the implicit scheme
% has much more numerical diffusion than the explicit scheme early in the
% simulation, but as the time increase, the difference is smaller. To
% verify this, we can plot the error or the breakthrough curves
%
n = size(pe,1);
pargs = {'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5]};
subplot(2,1,1),
   plot(1:n,e*100,'-o', pargs{:}),
   title('Percentage saturation discrepancy')
subplot(2,1,2),
   plot(1:n,pe(:,1),'-o',1:n,pi(:,1),'-s',pargs{:})
   legend('Explicit','Implicit','Location','NorthWest');
   title('Water breakthrough at heel'); axis tight

%%

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
