function [] = cleanupDialogue(precompDir)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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
