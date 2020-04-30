function data = computeDiagnostics(G, data, maxTOF, ix, precomp)
if nargin < 5
    precomp = [];
end
if nargin < 4 || isempty(ix)
    ix = 1:numel(data.states);
end

%  reverse engineer rock-struct so pore-volume matches ECLIPSE
if isfield(G.cells,'PORV')
    G.cells.volumes = G.cells.PORV;
    rock = struct('poro', ones(G.cells.num, 1));
else
    rock = data.rock;
end


vb = mrstVerbose && isempty(precomp);
dispif(vb, 'Computing diagnostics:     ');
for k = 1:numel(ix)
    if isempty(precomp)
        dispif(vb, '\b\b\b\b\b%3.0d %%', round(100*k/numel(ix)));
        st = data.states{ix(k)};
        W = data.wells{ix(k)};
        mrstVerbose('off');
        D = computeTOFandTracer(st, G, rock, 'wells', W, 'maxTOF', maxTOF, 'computeWellTOFs', true, ...
                                'processCycles', true, 'firstArrival', true);
        mrstVerbose(vb);
        % set tof to years
        data.diagnostics(ix(k)).D  = D;
        data.diagnostics(ix(k)).WP = computeWellPairs(st, G, rock, W, D );
    else
        data.diagnostics(k) = precomp{k}.diagnostics;
    end
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