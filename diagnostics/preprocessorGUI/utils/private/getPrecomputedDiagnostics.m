function precomp = getPrecomputedDiagnostics(casenm, steps, pdir)
% Function accepts either casenm for ECLIPSE input or problem structure for
% MRST input.


% if not given, assume precompted diagnositics lies in caseDir/mrst_diagnostics
if ischar(casenm) % ECLIPSE Input
    [caseDir, prefix] = fileparts(casenm);
elseif isstruct(casenm) % MRST Input
    caseDir = casenm.OutputHandlers.states.getDataPath;
    prefix = casenm.BaseName;
end

if nargin < 3
    pdir = fullfile(caseDir, 'mrst_diagnostics');
end

if ~isempty(dir(pdir))
    fn = @(n)fullfile(pdir, [prefix, sprintf('_diagn%0.4d.mat', n)]);
    fail = false;
    fprintf('Loading precomputed diagnostics %3.0d%%', 0)
    precomp = cell(numel(steps), 1);
    for k = 1:numel(steps)
        try
            precomp{k} = load(fn(steps(k)));
        catch
            fail = true;
            break
        end
        fprintf('\b\b\b\b%3.0d%%', round(100*k/numel(steps)))
    end
    fprintf(', done\n')
    if fail
        warning('Not able to load precomputed diagnostics')
        precomp = [];
    end
else
    precomp = [];
end
end
