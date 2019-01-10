function [] = precomputeDialogue(problem, precompDir)
sp = repmat(' ', [1, 8]);
    widestr = @(str)[sp, str, sp];
    
answ = questdlg(sprintf('Precompute diagnostics and write to:                                         \n   %s  \nThis might take up to 30 min',  precompDir), ...
             'Precompute for improved performance', ...
             widestr('OK'), widestr('skip'), widestr('OK'));
         
if strcmp(answ, widestr('OK'))
    if ischar(problem)
        % ECLIPSE
        processRestartDiagnostics(problem, 'outputdir', precompDir);
    elseif isstruct(problem)
        % If processing MRST 
        processStatesDiagnostics(problem, 'outputdir', precompDir);
    end
end
end