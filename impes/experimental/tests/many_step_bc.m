%%
assert ((exist('ROOTDIR', 'file') == 2) && isdir(ROOTDIR()), ...
        'MRST must be started prior to running this test.');

%%
clear

MODS = mrstModule;
VERB = mrstVerbose;
GRAV = gravity;

mrstModule  add   eclipse blackoil blackoiltransport
mrstModule('add', fileparts(mfilename('fullpath')));

mrstVerbose true
gravity     off

%%
ifluid = @(fn) ...
   initEclipseFluid(convertDeckUnits(readEclipseDeck( ...
   fullfile(ROOTDIR, 'params', 'fluid', 'Data', [fn, '.txt']))));

%%
cartDims = [100, 1];
physDims = [100, 1];

g    = computeGeometry(cartGrid(cartDims, physDims));
rock = struct('perm', repmat(0.1*darcy, [g.cells.num, 1]), ...
              'poro', repmat(1        , [g.cells.num, 1]));

%%
%fluid = ifluid('immiscibleoilgas');
fluid = ifluid('debug');

%%
bc = pside([], g, 'left' , 300*barsa, 'sat', [0, 1]); % Inject gas
bc = pside(bc, g, 'right', 100*barsa, 'sat', [1, 0]); % Produce oil

%%
T = computeTrans(g, rock);
x = initResSolComp(g, [], fluid, 100*barsa, [1, 0]); % Initially Oil-filled

%%
opts = {'bc', bc, 'verbose', true};

DT = 1*minute;
PV = poreVolume(g, rock);

x2 = x;
x3 = x;

X    = g.cells.centroids(:,1);
lbnd = min(g.faces.centroids(:,1));
ubnd = max(g.faces.centroids(:,1));

%%
for kk = 1:100,
   [x2, DT] = impesTPFA(x2, g, T, fluid, 10*DT, PV, opts{:}, ...
                 'ATol', 5.0e-11, 'RTol', 5.0e-13, 'urc', false);

   x3 = compTPFA(x3, g, rock, T, fluid, DT, opts{:});
   x3 = explicitTransport(x3, g, DT, rock, fluid, opts{:});

   figure(1)
   subplot(3, 1, 1),
   press = convertTo([x2.pressure, x3.pressure], barsa);
   plot(X, press(:,1), 'b-', X, press(:,2), 'g--');
   legend('impesTPFA', 'compTPFA', 'Location', 'Best');
   v = axis;
   axis([lbnd, ubnd, v(3:end)]);

   subplot(3, 1, 2),
   plot(X, x2.z(:,1), 'b-', X, x3.z(:,1), 'g--');
   legend('z_1 (IMPES)', 'z_1 (Seq Split)', 'Location', 'Best');
   v = axis;
   axis([lbnd, ubnd, v(3:end)]);

   subplot(3, 1, 3)
   plot(X, x2.z(:,2), 'b-', X, x3.z(:,2), 'g--');
   legend('z_2 (IMPES)', 'z_2 (Seq Split)', 'Location', 'Best');
   v = axis;
   axis([lbnd, ubnd, v(3:end)]);
end

%%
gravity(GRAV);  mrstVerbose(VERB);

mrstModule clear
mrstModule('add', MODS{:});
