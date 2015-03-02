function simtest

mrstModule add ad-fi

gravity on;

%% define parameters

celldim = [50, 1, 1];          % cells in Cartesian 3D grid
physdim = [5000, 3000, 100];   % physical dimensions of 3D grid
perm    = 100*milli*darcy;     % darcies
poro    = 0.1;                 % porosity
init_press  = 100 * barsa;     % Initial pressure (in barsa)

total_time  = 4 *  year;
nsteps      = 24; %12 * 20; % 48; %24; % one per month
dt = total_time / nsteps;

theta = (0/360) * 2 * pi;
slopedir = [1 0]; % sloping upwards (or downwards) along x-axis

%% Define grid and rock

[Gt, G] = topSurfaceGrid(computeGeometry(cartGrid(celldim, physdim)));

rock.perm = ones(G.cells.num, 1) * perm;
rock.poro = ones(G.cells.num, 1) * poro;

rock2D = averageRock(rock, Gt);


%% Define fluid properties

CO2 = CO2props;
temperature = 310; % degrees Kelvin; in Celsius it's 36.85

fluid.CO2.rho    = CO2.rho % @(p) CO2.rho    (p, temperature);
fluid.CO2.beta   = CO2.beta %@(p) CO2.beta   (p, temperature);
fluid.CO2.bder   = CO2.betader %@(p) CO2.betader(p, temperature);
fluid.CO2.mu     = 5.36108e-5; % in Pa/s at 310 deg. K
fluid.water.mu   = 6.5e-4;   % in Pa/s at 40 deg. C (313.15 K)
fluid.water.beta = 4.2e-10;
fluid.water.rho  = @(p, t) 1000 * (1 + fluid.water.beta * p); % first-order approximation used here

fluid_incomp = fluid;
fluid_incomp.water.rho = fluid.water.rho(init_press, temperature);
fluid_incomp.CO2.rho   = fluid.CO2.rho(init_press, temperature); %685.1531; % at 310 deg K, 100 barsa.

%fluid.muC = 7.858e-5; % in Pa/s at 40 deg. C (http://www.peacesoftware.de/einigewerte/co2_e.html)
% fluid.rhoC = 839;   % at 40 deg. C and 20 MPa                                                                     

%% load VE correction polynomials (only used when verifying mass balance during simulation below)
VEpoly = VEpolys;

%% define wells

%computing a reasonable rate:
if false
    cell_vol = poro * (prod(physdim) / prod(celldim)); % volume in m^3
    rate = (1/20) * cell_vol * meter^3 / day;    % should take 20 days to inject
                                                 % volume equivalent to a cell.
    mass_rate = rate * fluid.rhoC(init_press);
else
    mass_rate = 1e6 * kilo * kilogram / year; % in tonnes - chosen to be
                                              % comparable to Utsira
    % mass_rate = mass_rate * 50; % HACK @@
    rate = mass_rate / fluid.CO2.rho(init_press, temperature);
end

%mass_rate = 0; rate = 0; % @@ hack

% defining an injector right in the middle
W_incomp = verticalWell([], G, rock, (celldim(1)+0)/2, 1, 1:celldim(3), ...
                        'Type', 'Rate', 'Val', rate, 'Radius', 0.1, 'Name', 'I', ...
                        'Comp_i', [1 0]);
Wt_incomp = convertwellsVE(W_incomp, G, Gt, rock2D);

W = verticalWell([], G, rock, (celldim(1)+0)/2, 1, 1:celldim(3), ...
                 'Type', 'Rate', 'Val', mass_rate, 'Radius', 0.1, 'Name', 'I', ...
                 'Comp_i', [1 0]);

Wt = convertwellsVE(W, G, Gt, rock2D);

%% Define helper structure

% computes: - face transmissibilities (all, and internal)
%           - pore volumes
%           - divergence matrix C
%           - function for computing internal face averages of some field
%           - function for projecting cell values to upstream faces
%           - neighbor relations (N)

% calling `setupSimComp` instead of `setupSimCompVe`, since there's at least
% one bug in the latter, and since it's apparently prepared for dealing with
% saturations/pore volume rather than porosity.  Instead, I call the regular
% `setupSimComp function, but modify the permeability of the rock first, in
% order to make it not the vertical average, but the vertical sum.
rock2D.perm = rock2D.perm .* Gt.cells.H; 
s = setupSimComp(Gt, rock2D);

s_incomp = s;



%% Set initial state

xmin = min(Gt.faces.centroids(:,1));
dheight = (Gt.cells.centroids(:,1) - xmin) * sin(theta); %@ assumes sloping direction in x here!

% hydrostatic pressure
x0.temperature = temperature; % not expected to change.
x0.pressure = init_press*ones(Gt.cells.num, 1) - fluid.water.rho(init_press) * norm(gravity) * dheight;

% x0.pressure = init_press * ones(Gt.cells.num, 1); % equal initial pressure
%                                               % everywhere
x0.h = zeros(Gt.cells.num, 1);                % initially filled with water everywhere

%x0.h(10:15) = 50;  %% @@ HACK
%x0.h(20:25) = 50;
% x0.h(40:45) = 50;  %% @@ HACK

x0_incomp = x0;


%% Define boundary conditions

%bc_val = init_press;
% bc = pside([] ,Gt, 'LEFT' , bc_val);
% bc = pside(bc, Gt, 'RIGHT', bc_val);

bc = pside([], Gt, 'LEFT', x0.pressure(1));
bc = pside(bc, Gt, 'RIGHT', x0.pressure(end));

% bc = fluxside([], Gt, 'LEFT',  0); no flow
% bc = fluxside(bc, Gt, 'RIGHT', 0); %no flow

%% Define system

system.stepFunction = @(x0, x, meta, dt, W, G, system) stepCO2(x0, x, meta, ...
                                                  dt, W, G, system, fluid);
system.getEquations = @(state0, state, dt, G, W, s, f) ...
    eqsfiCO2compressibleSimplest(state0, state, dt, G, W, s, f, bc, theta, slopedir);
system.s = s;                            % storing helper structure in system
system.nonlinear.maxIterations = 25;     % queried inside 'stepCO2'
system.nonlinear.changeWells = false;    % queried in 'solvefiADI'
system.nonlinear.bhpcontrols = false;    % queried in 'solvefiADI'
system.nonlinear.relaxation  = 1;        % queried in 'solvefiADI' but not used!
system.nonlinear.tolMB = 1e-7;           % queried in 'getCO2Convergence'
system.nonlinear.tolCNV = 1e-3;          % queried in 'getCO2Convergence'

system_incomp = system;
system_incomp.stepFunction = @(x0, x, meta, dt, W, G, system) stepCO2(x0, x, meta, dt, W, G, system, fluid_incomp);
system_incomp.getEquations = @(state0, state, dt, G, W, s, f) ...
    eqsfiCO2simplest(state0, state, dt, G, W, s, f, bc, theta, slopedir);

%% Solve system

x = x0;
x_incomp = x0;

x_incomp_escaped = 0;
x_escaped = 0;

xc=Gt.cells.centroids(:,1);
figure(2);
for tstep = 1:nsteps
    [x, its]        = solvefiADI(x,        dt, Wt, Gt, system);
    [x_incomp, its] = solvefiADI(x_incomp, dt, Wt_incomp, Gt, system_incomp);

    subplot(5,1,1),cla, hold on
      dualplot(xc/1e3, x.pressure/1e6, x_incomp.pressure/1e6, 'km', 'MPa');

    subplot(5,1,2), cla, hold on
      ipress_incomp = ...
          x_incomp.pressure ...
          + fluid_incomp.CO2.rho * norm(gravity) * cos(theta) ...
          * x_incomp.h;
      ipress = ...
          x.pressure ...
          + fluid.CO2.rho(x.pressure, x.temperature) * norm(gravity) * cos(theta) .* ...
          VEpoly.intAlpha(x.h, fluid.CO2.rho(x.pressure, x.temperature), fluid.CO2.beta(x.pressure, x.temperature), ...
                     fluid.CO2.bder(x.pressure, x.temperature), theta);
      dualplot(xc/1e3, ipress/1e6, ipress_incomp/1e6, 'km','MPa');
    subplot(5,1,3), cla, hold on
      bpress_incomp = ...
          ipress_incomp + fluid_incomp.water.rho * norm(gravity) * cos(theta) * ...
          (Gt.cells.H - x_incomp.h);
      bpress = ...
          ipress + fluid_incomp.water.rho * norm(gravity) * cos(theta) * ...
          (Gt.cells.H - x.h);
      dualplot(xc/1e3, bpress/1e6, bpress_incomp/1e6, 'km', 'MPa');
    subplot(5,1,4),cla, hold on
      dualplot(xc/1e3, x.h, x_incomp.h, 'km', 'meters');
      % plot(xc,x.h, 'r');     
      % plot(xc, x_incomp.h, 'b');
      axis([0 physdim(1)/1e3 0 physdim(3)]); set(gca, 'YDir', 'reverse');

    subplot(5,1,5), cla, hold on
      dualplot(xc, fluid.CO2.rho(x.pressure, x.temperature), fluid_incomp.CO2.rho * ones(size(x.pressure, x.temperature)), ...
               'km', 'kg/m3');
      % plot(xc, fluid.rhoC(x.pressure), 'r');
      % plot(xc, fluid_incomp.rhoC * ones(size(x.pressure)), 'b');

    x_incomp_escaped = x_incomp_escaped + sum(x_incomp.info.outflowCO2) * dt;
    x_escaped        = x_escaped        + sum(x.info.outflowCO2) * dt;   
    
    
    rhoC = fluid.CO2.rho(x.pressure, x.temperature);
    VErhoC = rhoC .* VEpoly.intAlpha(x.h, rhoC, fluid.CO2.beta(x.pressure,x.temperature), ...
                                fluid.CO2.bder(x.pressure,x.temperature), theta);
    co2_in_grid        = sum(Gt.cells.volumes .* rock.poro .* VErhoC);
    co2_in_grid_incomp = x_incomp.h' * Gt.cells.volumes * poro * fluid_incomp.CO2.rho;
    
    fprintf('\n---------\n\n');
    fprintf('Injected (incomp) so far: %e \n', tstep * dt * Wt_incomp.val * fluid_incomp.CO2.rho);
    fprintf('Contained in grid: %e\n', co2_in_grid_incomp);
    fprintf('Escaped: %e\n', x_incomp_escaped.val * fluid_incomp.CO2.rho);
    fprintf('Balance: %e\n', x_incomp_escaped.val * fluid_incomp.CO2.rho + ...
            co2_in_grid_incomp - tstep*dt*Wt_incomp.val * fluid_incomp.CO2.rho);
    
    fprintf('\n');
    fprintf('Injected (comp) so far: %e\n', tstep * dt * Wt.val);
    fprintf('Contained in grid: %e\n', co2_in_grid);
    fprintf('Escaped: %e\n', x_escaped.val);
    fprintf('Balance: %e\n', x_escaped.val + co2_in_grid - tstep*dt*Wt.val);

    %plotCellData(Gt, x.h); % colorbar('horiz', caxis([0 1]));
    drawnow;
    
    %print('-dpng', '-r200', sprintf('year %03d.png', tstep*dt/year));
    %pause;
    %pause(0.05); %pause(0.05);
end

end % end of main function
    
function dualplot(x, comp, incomp, xl, yl)
    plot(x, comp, 'r');
    plot(x, incomp, 'b');
    xlabel(xl);
    ylabel(yl);
end

%% Hypotheses
% * Zero incliniation
% * No residual saturations
% * Flat top and bottom
% * incompressible
% * zero entry pressure
% * sharp interface
% * constant permeability and porosity
% * constant pressure boundary conditions
