function debug_gravity()
%% Define geometry

gravity([0, 0, 10]), gravity on

cartDims = [1, 1, 10];
physDims = [1,1,200]*meter;

G = computeGeometry(cartGrid(cartDims, physDims));


%% Define medium and fluid properties

rock = struct('perm', repmat(1, [G.cells.num, 1])*darcy, ...
              'poro', repmat(1, [G.cells.num, 1]));

%%
fluid = initSimpleCompFluid('mu', [1, 1, 1], 'rho_ref', [1000, 1000, 600], ...
                            'press_ref', 100*barsa, 'n', [1, 1, 1], ...
                            'c', [1, 9, 9]/(100*barsa));

%% Initial reservoir state and driving forces

initP     = 100*barsa;
initcomp  = [1, 0];

% Boundary conditions -----------------------------------------------------
%
injcomp = [0, 1];

bc = [];
%bc = addBC(bc, 1, 'pressure', 100*barsa, 'sat', injcomp);
deck.SOLUTION.EQUIL=[100   100*barsa    300    0     100     0  0 0 0 0 0];
deck.PROPS = [];
deck.REGIONS = [];
%state1 = initResSolComp(G, [], fluid, initP, initcomp);
state1 = initEclipseState(G,deck,fluid);
%state1.pressure=state1.pressure+1*barsa;
state2 = state1;



%% Solvers

htrans = computeTrans(G, rock);
trans  = 1 ./ accumarray(G.cells.faces(:,1), ...
                         1 ./ htrans, [G.faces.num, 1]);

%%{
pvol    = poreVolume(G, rock);
psolve1 = @(state, dt, vd) ...
   impesTPFA(state, G, trans, fluid, dt, pvol, ...
             'bc', bc, 'RTol', 5.0e-8,'EstimateTimeStep',true,'UpdateMass',true,'DynamicMobility',false,'face_update',false);

%{
% Cheating.  Solve an incompressible problem...
fluid2 = initSimpleFluid('mu', [1, 1], 'rho', [700, 1000], 'n', [2, 2]);
psolve2 = @(state, dt, varargin) ...
   incompTPFA(state, G, htrans, fluid2, 'bc', bc);

tsolve2 = @(state, dt) ...
   explicitTransport(state, G, dt, rock, fluid2, 'bc', bc);
%}
%%{

psolve2 = @(state, dt, vd) ...
   compTPFA(state, G, rock, htrans, fluid, dt, ...
            'bc', bc, 'volume_correction', vd);

tsolve2 = @(state, dt) expl_bo_transport(state, G, dt, rock, fluid, bc);
%}

p_ana = @(x) -log(-comp*barsa*x + 1) / (comp * barsa);

%% Time loop
t  = 0;
dt = 0.1*year;
X  = G.cells.centroids(:,3);

vol_d = @(u, dt) (sum(u,2) - 1) .* pvol ./ dt;
figure(1);figure(2)
for k = 1:200,
   [state1, dt] = psolve1(state1, 2*dt, zeros([G.cells.num, 1]));
%{
   [u, u, u, u] = fluid.pvt(state2.pressure, state2.z);                %#ok
   state2       = psolve2(state2,   dt, vol_d(u, dt));
   state2       = tsolve2(state2,   dt);
%}
   t = t + dt;
   disp(['Elapsed time', num2str(t/year),' year'])
   % Plotting below -------------------------------------------------------

   set(0, 'CurrentFigure', 1)
   subplot(3,1,1),
   plot(X, convertTo([state1.pressure , state2.pressure], barsa), '.')
   title('pressure')
   legend('impes','compTPFA exp')
   subplot(3,1,2),
   plot(X, convertTo( state1.pressure - state2.pressure , barsa), '.')
   title('Pressure difference')
   subplot(3,1,3),
   plot(X, state1.z - state2.z , '.')
   title('dz difference')

   set(0, 'CurrentFigure', 2)
   %plot(X, state1.z, '.', X, state2.z, '-x')
   plot(X, state1.z, '*')
   hold on;
   plot(X, state2.z, '-x')
   hold off
   title('z')
   %legend('impes z_1','impes z_2','compTPFA exp')
   pause(0.25)
end

%--------------------------------------------------------------------------

function state = expl_bo_transport(state, G, dt, rock, fluid, bc)
%MODS = mrstModule;
%mrstModule add blackoiltransport

state = implicitTransport(state, G, dt, rock, fluid, 'bc', bc);

%mrstModule('reset', MODS{:});

