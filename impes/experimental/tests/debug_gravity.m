function debug_gravity()
%% Define geometry

gravity([1, 0, 0]), gravity on

cartDims = 200;
physDims = 200*meter;

G = cartGrid1D(cartDims, physDims);


%% Define medium and fluid properties

rock = struct('perm', repmat(1, [G.cells.num, 1]), ...
              'poro', repmat(1, [G.cells.num, 1]));

%%
fluid = initSimpleCompFluid('mu', [1, 1], 'rho_ref', [700, 1000], ...
                            'press_ref', 100*barsa, 'n', [1, 1], ...
                            'c', 1*[9e-5, 9e-6]/barsa);

%% Initial reservoir state and driving forces

initP     = 100*barsa;
initcomp  = [1, 0];

% Boundary conditions -----------------------------------------------------
%
injcomp = [0, 1];

bc = [];
bc = addBC(bc, 1, 'pressure', 100*barsa, 'sat', injcomp);

state1 = initResSolComp(G, [], fluid, initP, initcomp);
state2 = state1;


%% Solvers

htrans = computeTrans(G, rock);

%%{
pvol    = poreVolume(G, rock);
psolve1 = @(state, dt, vd) ...
   impesTPFA(state, G, htrans, fluid, dt, pvol, ...
             'bc', bc, 'RTol', 5.0e-12);

%{
% Cheating.  Solve an incompressible problem...
fluid2 = initSimpleFluid('mu', [1, 1], 'rho', [700, 1000], 'n', [2, 2]);
psolve2 = @(state, dt, varargin) ...
   incompTPFA(state, G, htrans, fluid2, 'bc', bc);

tsolve2 = @(state, dt) ...
   explicitTransport(state, G, dt, rock, fluid2, 'bc', bc);
%}
%{

psolve2 = @(state, dt, vd) ...
   compTPFA(state, G, rock, htrans, fluid, dt, ...
            'bc', bc, 'volume_correction', vd);

tsolve2 = @(state, dt) expl_bo_transport(state, G, dt, rock, fluid, bc);
%}

p_ana = @(x) -log(-comp*barsa*x + 1) / (comp * barsa);

%% Time loop
t  = 0;
dt = 1*day;
X  = G.cells.centroids;

vol_d = @(u, dt) (sum(u,2) - 1) .* pvol ./ dt;

for k = 1:200,
   [state1, dt] = psolve1(state1, 2*dt, zeros([G.cells.num, 1]));

   [u, u, u, u] = fluid.pvt(state2.pressure, state2.z);                %#ok
   %state2       = psolve2(state2,   dt, vol_d(u, dt));
   %state2       = tsolve2(state2,   dt);

   t = t + dt;

   % Plotting below -------------------------------------------------------

   set(0, 'CurrentFigure', 1)
   subplot(2,1,1),
   plot(X, convertTo([state1.pressure , state2.pressure], barsa), '.')
   subplot(2,1,2),
   plot(X, convertTo( state1.pressure - state2.pressure , barsa), '.')

   set(0, 'CurrentFigure', 2)
   plot(X, state1.z, '.', X, state2.z, '-x')

   pause(0.25)
end

%--------------------------------------------------------------------------

function state = expl_bo_transport(state, G, dt, rock, fluid, bc)
MODS = mrstModule;
mrstModule add blackoiltransport

state = implicitTransport(state, G, dt, rock, fluid, 'bc', bc);

mrstModule('reset', MODS{:});

