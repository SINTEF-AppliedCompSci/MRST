% interactiveDiagnosticsFromEclipseOutput
mrstModule add mrst-gui ad-core diagnostics

prefix =  '/Users/steink/data/opm-data/spe9/eclipse-simulation/SPE9_CP';

% Read INIT/EGRID-files and construct MRST-grid
init = readEclipseOutputFileUnFmt([prefix, '.INIT']);
egrid = readEclipseOutputFileUnFmt([prefix, '.EGRID']);
[G, rock] = eclOut2mrst(init, egrid);
G = computeGeometry(G);

% Convert restart step 10 to mrst-state:
states = convertRestartToStates(prefix, G, 'steps', 10);
state = states{1};

% Well struct from wellSol
W = state.wellSol;

% If flux field is not given in Eclipse restart, state will not contain
% fluxes, and interactiveDiagnostics will solve a pressure equation. We
% adjust W such that it the control settings are compatible with the
% incompressible solvers, hence we set W(1) to be controlled by bhp, and
% W(2:end) to be controlled by rate. We set the rate values to the
% reservoir volume rates computed by ECLIPSE:
[W(1).type, W(1).val] = deal('bhp', W(1).bhp);
for wnum = 2:numel(W)
    [W(wnum).type, W(wnum).val] = deal('rate', W(wnum).resv);
end

% run Interactive diagnostics
interactiveDiagnostics(G, rock, W, 'state', state);
figure(1), daspect([1 1 .3]), view([1 -.5 1])