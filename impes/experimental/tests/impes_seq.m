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
cartDims = [40, 1];
physDims = [5, 1];

g    = computeGeometry(cartGrid(cartDims, physDims));
rock = struct('perm', repmat(0.1*darcy, [g.cells.num, 1]), ...
              'poro', repmat(1        , [g.cells.num, 1]));

%%
fluid = ifluid('immiscibleoilgas');
% fluid = ifluid('debug');

%%
W = [];
%{
W = addWell(W, g, rock, 1, 'Type', 'bhp', 'Val', 300*barsa, ...
            'Radius', 10*centi*meter, 'Dir', 'z',           ...
            'InnerProduct', 'ip_tpf', 'Comp_i', [0, 1]); % Inject gas
%}
%%{
W = addWell(W, g, rock, 1, 'Type', 'rate', 'Val', 1*meter^3/day, ...
            'Radius', 10*centi*meter, 'Dir', 'z',                ...
            'InnerProduct', 'ip_tpf', 'Comp_i', [0, 1]); % Inject gas
%}
W = addWell(W, g, rock, g.cells.num, 'Type', 'bhp', 'Val', 100*barsa, ...
            'Radius', 10*centi*meter, 'Dir', 'z',                     ...
            'InnerProduct', 'ip_tpf', 'Comp_i', [1, 0]); % Produce oil

%%
T = computeTrans(g, rock);
x = initResSolComp(g, W, fluid, 100*barsa, [1, 0]); % Initially Oil-filled

trans = 1 ./ accumarray(g.cells.faces(:,1), 1 ./ T, [g.faces.num, 1]);

%%
opts = {'wells', W, 'verbose', true};

DT = 0.5*minute; %3*hour; %10*minute;
PV = poreVolume(g, rock);

x2 = x;
x3 = x;

X    = g.cells.centroids(:,1);
lbnd = min(g.faces.centroids(:,1));
ubnd = max(g.faces.centroids(:,1));

%%
for kk = 1:20,%300,
   [x2, DT] = impesTPFA(x2, g, trans, fluid, 2*DT, PV, opts{:}, ...
                        'UpdateMass', false, 'ATol', 5.0e-11, ...
                        'RTol', 5.0e-13);

   x2 = explicitTransport(x2, g, DT, rock, fluid, opts{:});

%%{
   x3 = compTPFA(x3, g, rock, T, fluid, DT, opts{:});
   x3 = explicitTransport(x3, g, DT, rock, fluid, opts{:});
%}

   figure(1)
   subplot(2, 1, 1),
   plot(X, x2.z(:,1), 'b-', X, x3.z(:,1), 'g--');
   legend('z_1 (Combined)', 'z_1 (Seq Split)', 'Location', 'Best');
   v = axis;
   axis([lbnd, ubnd, v(3:end)]);

   subplot(2, 1, 2)
   plot(X, x2.z(:,2), 'b-', X, x3.z(:,2), 'g--');
   legend('z_2 (Combined)', 'z_2 (Seq Split)', 'Location', 'Best');
   v = axis;
   axis([lbnd, ubnd, v(3:end)]);

   pause(0.2)
end

%%
gravity(GRAV);  mrstVerbose(VERB);

mrstModule clear
mrstModule('add', MODS{:});
