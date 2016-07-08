%% Read/convert ECLIPSE output and visualize using interactive diagnostics
% This example reads ECLIPSE unformatted output files from a simulation
% based on the SPE-9 benchmark. The data is converted to MRST grid and
% soluion structures and used as input to the interactive diagnostics
% functionality.

mrstModule add mrst-gui ad-core diagnostics

if ~ makeSPE9OutputAvailable,
   error('SPE9Download:Failure', ...
         'Failed to download ECLIPSE output for SPE-9 benchmark case');
end

prefix = fullfile(getDatasetPath('spe9'), 'Simulation-Output', 'SPE9_CP');

%% Read INIT/EGRID-files and construct MRST-grid
init = readEclipseOutputFileUnFmt([prefix, '.INIT']);
egrid = readEclipseOutputFileUnFmt([prefix, '.EGRID']);
[G, rock] = eclOut2mrst(init, egrid);
G = computeGeometry(G);

%% Read restart step 10 and convert to mrst-state
step = 10;
states = convertRestartToStates(prefix, G, 'steps', 10);
state = states{1};

%% Set up a compatible well structure (W)
% The wellSol-field of state contains sufficient data to work as a
% well-struct (W).
W = state.wellSol;

%%
% If flux field is not given in ECLIPSE restart (as here), the converted
% state will not contain fluxes either. In this case,
% interactiveDiagnostics will need to solve an incompressible pressure
% equation in order to reconstruct a reasonable flux field. We therefore
% need to adjust |W| such that its control settings are compatible with the
% incompressible solvers. Here we set |W(1)| (the injector) to be
% controlled by bottom-hole pressure target (BHP), and |W(2:end)| (the
% producers) to be controlled by rate. We set the rate values to the
% reservoir volume rates computed by ECLIPSE (resv).

[W(1).type, W(1).val] = deal('bhp', W(1).bhp);
for wnum = 2:numel(W)
    [W(wnum).type, W(wnum).val] = deal('rate', W(wnum).resv);
end

%% Run Interactive diagnostics
close all, interactiveDiagnostics(G, rock, W, 'state', state);
figure(1), daspect([1, 1, 0.3]), view([1, -0.5, 1])
