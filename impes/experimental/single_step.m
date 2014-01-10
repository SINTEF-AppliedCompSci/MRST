%%
assert ((exist('ROOTDIR', 'file') == 2) && isdir(ROOTDIR()), ...
        'MRST must be started prior to running this test.');

%%
clear

MODS = mrstModule;
VERB = mrstVerbose;
GRAV = gravity;

mrstModule  add   eclipse blackoil blackoiltransport
% mrstModule('add', fileparts(mfilename('fullpath')));

mrstVerbose true
gravity     off

%%
ifluid = @(fn) ...
   initEclipseFluid(convertDeckUnits(readEclipseDeck( ...
   fullfile(ROOTDIR, 'params', 'fluid', 'Data', [fn, '.txt']))));

%%
cartDims = [20, 1];
physDims = cartDims;

g    = computeGeometry(cartGrid(cartDims, physDims));
rock = struct('perm', repmat(0.1*darcy, [g.cells.num, 1]), ...
              'poro', repmat(1        , [g.cells.num, 1]));

%%
fluid = ifluid('immiscibleoilgas');

%%
bc = pside([], g, 'left' , 300*barsa, 'sat', [0, 1]); % Inject gas
bc = pside(bc, g, 'right', 100*barsa, 'sat', [1, 0]); % Produce oil

%%
T = computeTrans(g, rock);
x = initResSolComp(g, [], fluid, 100*barsa, [1, 0]); % Initially Oil-filled

%%
opts = {'bc', bc, 'verbose', true};

DT = 1*minute;

%%
x2 = impesTPFA(x, g, T, fluid, DT, poreVolume(g, rock), opts{:}, ...
               'ATol', 5.0e-11, 'RTol', 5.0e-10);

%%
x3 = compTPFA(x, g, rock, T, fluid, DT, opts{:});
x3 = explicitTransport(x3, g, DT, rock, fluid, opts{:});

%%
gravity(GRAV);  mrstVerbose(VERB);

mrstModule clear
mrstModule('add', MODS{:});

%%
X    = g.cells.centroids(:,1);
lbnd = min(g.faces.centroids(:,1));
ubnd = max(g.faces.centroids(:,1));

figure(1)
subplot(3, 1, 1),
press = convertTo([x2.pressure, x3.pressure], barsa);
plot(X, press(:,1), 'g-', X, press(:,2), 'r--');
legend('impesTPFA', 'compTPFA', 'Location', 'Best');
v = axis;
axis([lbnd, ubnd, v(3:end)]);

subplot(3, 1, 2),
plot(X, x2.z(:,1), 'g-', X, x3.z(:,1), 'r--');
legend('z_1 (IMPES)', 'z_1 (Seq Split)', 'Location', 'Best');
v = axis;
axis([lbnd, ubnd, v(3:end)]);

subplot(3, 1, 3)
plot(X, x2.z(:,2), 'g-', X, x3.z(:,2), 'r--');
legend('z_2 (IMPES)', 'z_2 (Seq Split)', 'Location', 'Best');
v = axis;
axis([lbnd, ubnd, v(3:end)]);
