function [] = cleanupDialogue(precompDir)
[~,pd] = fileparts(precompDir);
assert(strcmp(pd, 'mrst_diagnostics'), 'Attempt to run cleanup on unexpected directory')
if ~isempty(ls(precompDir))
    sp = repmat(' ', [1, 8]);
    widestr = @(str)[sp, str, sp];
    dd = dir(precompDir);
    nf = nnz([dd.bytes]);
    answ = widestr('OK');
    if nf > 0 
        answ = questdlg(sprintf('Directory: \n %s \n contains %d non-empty files/folders. \n Really remove?', precompDir, nf), ...
                        'Remove pre-computed diagnostics', ...
                        widestr('OK'), widestr('cancel'), widestr('OK'));
    end
    if strcmp(answ, widestr('OK'))
        rmdir(precompDir, 's');
    end
end
end