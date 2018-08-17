function data = computeDiagnostics(G, data, maxTOF, ix)
if nargin < 4 || isempty(ix)
    ix = 1:numel(data.states);
end

%  reverse engineer rock-struct so pore-volume matches ECLIPSE
G.cells.volumes = G.PORV;
rock = struct('poro', ones(G.cells.num, 1));

vb = mrstVerbose;
dispif(vb, 'Computing diagnostics:     ');
for k = 1:numel(ix)
    
    dispif(vb, '\b\b\b\b\b%3.0d %%', round(100*k/numel(ix)));
    st = data.states{ix(k)};
    mrstVerbose('off');
    D = computeTOFandTracerFirstArrival(st, G, rock, 'wells', st.wellSol, 'maxTOF', maxTOF, 'computeWellTOFs', true, ...
        'processCycles', true);
    mrstVerbose(vb);
    % set tof to years
    data.diagnostics(ix(k)).D  = D;
    data.diagnostics(ix(k)).WP = computeWellPairs(st, G, rock, st.wellSol, D );
end
dispif(vb, ', done\n');

% compute well-communication matrix - average over time-steps
dt = data.time.cur - data.time.prev;
com = 0;
for k = 1:numel(dt)
    salloc = cellfun(@sum, {data.diagnostics(k).WP.inj.alloc}, 'UniformOutput',false);
    salloc = vertcat(salloc{:});
    com = com + dt(k)*salloc/sum(dt); 
end
 
data.wellComunication = abs(com);% > 0.01*totalloc/n;
end