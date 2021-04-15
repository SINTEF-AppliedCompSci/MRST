%% Pressure Solver with capillary pressure:
% Here, we demonstrate the effect of capillary pressure on the solution of
% a two-phase oil-water problem. We solve the two-phase pressure equation
%
% $$\nabla\cdot v = q, \qquad v=\textbf{--}\lambda_t K\nabla p,$$
%
% where v is the Darcy velocity (total velocity) and lambda_t is the
% total mobility, which depends on the water saturation S.
%
% The saturation equation (conservation of the water phase) is given as:
%
% $$ \phi \frac{\partial S}{\partial t} +
%     \nabla \cdot (f_w(S)(v + K\lambda_o \nabla p_c)) = q_w$$
%
% This tutorial shows a 2D case with homogeneous permeability and
% porosity and linear capillary pressure curve and is based on the
% <matlab:edit('incompTutorialWells2ph.m') incompTutorialWells2ph>  example.
%
mrstModule add incomp
verbose = true;

%% Construct simple Cartesian test case
nx = 40; ny = 40; nz = 1;
G         = cartGrid([nx ny nz]);
G         = computeGeometry(G);
rock      = makeRock(G, 100*milli*darcy, 0.3);
hT  = computeTrans(G, rock, 'Verbose', verbose);

%% Define fluid and capillary pressure curve
% We define the relative permeability and the capillary pressure in form of
% tables, and let the relative permeability curves be quadratic and the
% capillary function linear. The strength of the capillary pressure is
% decided by cap_scale. The capillary pressure is defined in the
% non-wetting phase, i.e. $$ p_c = p_{nw} - p_w $$.
pc_form = 'nonwetting';
cap_scale = 10;
x = linspace(0, 1, 11) .';
y = linspace(1, 0, 11) .';
[kr, pc]  = tabulatedSatFunc([x, x.^2, y.^2, y.*cap_scale*barsa]);

clf,
subplot(1,2,1), plot(x,kr(x),'LineWidth',2); title('Relative permeability');
subplot(1,2,2), plot(x,pc(x),'LineWidth',2); title('Capillary function');

% Define constant properties for viscosity and density
props = constantProperties([   1,  10] .* centi*poise, ...
                           [1000, 700] .* kilogram/meter^3);

%%
% Here we put together a valid fluid object from the above defined
% functions. To read more about the fluid structure write
% help fluid_structure in MRST. First make a fluid without capillary
% pressure
fluid = struct('properties', props                  , ...
               'saturation', @(x, varargin)    x.s  , ...
               'relperm'   , kr);

%%
% Then make another fluid object identical to the one above except for the
% capillary pressure term 'pc'.
fluid_pc = struct('properties', props                  , ...
                  'saturation', @(x, varargin)    x.s  , ...
                  'relperm'   , kr                     , ...
                  'pc'        , @(x, varargin) pc(x.s));

%% Plot the pc-curve
% Make a dummy state/solution structure to plot the pc curve since
% 'fluid.pc' demands state as an input
xDummy   = initState(G, [], [0, 1]);
xDummy.s = linspace(0, 1, numel(xDummy.s))'; ...
pc = convertTo(fluid_pc.pc(xDummy), barsa);

clf
plot(xDummy.s, pc, 'LineWidth',2);
xlabel('s_w'); ylabel('pc [bar]');
title('Capillary pressure curve')


%% Set wells and ininitialize reservoir state
rate = 0.5*meter^3/day;
bhp  = 1*barsa;
W = verticalWell([], G, rock, 1, 1, 1:nz,          ...
                 'InnerProduct', 'ip_tpf', ...
                 'Type', 'rate', 'Val', rate, ...
                 'Radius', .1, 'Name', 'I', 'Comp_i', [1 0]);
W = verticalWell(W, G, rock, nx, ny, 1:nz,     ...
                 'InnerProduct', 'ip_tpf', ...
                 'Type','bhp', 'Val', bhp, ...
                 'Radius', .1, 'Dir', 'x', 'Name', 'P', 'Comp_i', [0 1]);

rSol    = initState(G, W, 0, [0.2, 0.8]);
rSol_pc = initState(G, W, 0, [0.2, 0.8]);
gravity off
verbose = false;

%% Set up pressure and transport solvers
% This example uses an implicit transport solver, an explicit solver can be
% used if the time step restriction for the parabolic term is less than for
% the hyperbolic term. This is the case if 'cap_scale' is small. We let
% 'fluid' be a parameter in 'psolve' and 'tsolve' so that we can use the
% solvers for simulation both with and without capillary pressure by
% supplying different fluid objects. For this case we use the verbose =
% false for the transport solver. If more information about the convergence
% of the method is required; use verbose = true.

psolve  = @(state, fluid) ...
    incompTPFA(state, G, hT, fluid, 'wells', W);
tsolve  = @(state, dT, fluid) ...
    implicitTransport(state, G, dT, rock, fluid, 'wells', W, 'verbose', verbose);
%
% Alternatively we could have defined an explicit transport solver by
% tsolve = @(state, dT) ...
%     explicitTransport(state, G, dT, rock, fluid, 'wells', W, 'verbose', verbose);

%% Solve initial pressure in reservoir
% Observe that we supply different fluid objects for the two solutions, one
% with capillary pressure and one without.
rSol    = psolve(rSol, fluid);
rSol_pc = psolve(rSol_pc, fluid_pc);

%% Transport loop
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation.
T      = 300*day();
dT     = T/15;
dTplot = 100*day();  % plot only every 100th day
N      = fix(T/dTplot);
pv     = poreVolume(G,rock);

%% Start the main loop
t  = 0; plotNo = 1;
h1 = 'No pc - '; h2 = 'Linear pc - ';
e = []; p_org = []; p_pc = [];
clf;

while t < T,
   % TRANSPORT SOLVE
   rSol    = tsolve(rSol, dT, fluid);
   rSol_pc = tsolve(rSol_pc, dT, fluid_pc);

   % Check for inconsistent saturations
   s = [rSol.s(:,1); rSol_pc.s(:,1)];
   assert(max(s) < 1+eps && min(s) > -eps);

   % Update solution of pressure equation.
   rSol    = psolve(rSol,    fluid);
   rSol_pc = psolve(rSol_pc, fluid_pc);

   % Measure water saturation in production cells in saturation
   e = [e; sum(abs(rSol.s(:,1) - rSol_pc.s(:,1)).*pv)/sum(pv)]; %#ok
   p_org = [p_org; rSol.s(W(2).cells,1)' ];                 %#ok
   p_pc  = [p_pc; rSol_pc.s(W(2).cells,1)'];                 %#ok

   % Increase time and continue if we do not want to plot saturations
   t = t + dT;
   if ( t < plotNo*dTplot && t <T), continue, end

   % Plot saturation
   heading = [num2str(convertTo(t,day)),  ' days'];
   r = 0.01;
   subplot('position',[(plotNo-1)/N+r, 0.50, 1/N-2*r, 0.48]), cla
   plotCellData(G, rSol.s(:,1),'EdgeColor','none');
   caxis([0 1]), view(60,50), axis tight off, set(gca,'dataasp',[12 12 1]),
   title([h1 heading])

   subplot('position',[(plotNo-1)/N+r, 0.02, 1/N-2*r, 0.48]), cla
   plotCellData(G, rSol_pc.s(:,1),'EdgeColor','none');
   caxis([0 1]), view(60,50), axis tight off, set(gca,'dataasp',[12 12 1]),
   title([h2 heading])
   drawnow

   plotNo = plotNo+1;

end

%% Plot water breakthrough at heel
% As we clearly see from the plots in the figure, the simulation with
% capillary pressure has much more diffusion than the simulation without
% capillary pressure. This is confirmed by the water breakthrough curve.
%
clf
n = numel(p_org(:,1));
plot(1:n,p_org(:,1),'-o',1:n,p_pc(:,1),'--s','MarkerSize',8,'MarkerFaceColor',[.7 .7 .7])
legend('No capillary pressure','Linear capillary pressure','Location','Best');
axis tight
title('Water breakthrough at heel');

%% Copyright notice

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
