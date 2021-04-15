function [] = precomputeDialogue(casenm, precompDir)
%Undocumentd Utility Function

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

sp = repmat(' ', [1, 8]);
    widestr = @(str)[sp, str, sp];
    
answ = questdlg(sprintf('Precompute diagnostics and write to:                                         \n   %s  \nThis might take up to 30 min',  precompDir), ...
             'Precompute for improved performance', ...
             widestr('OK'), widestr('skip'), widestr('OK'));
         
if strcmp(answ, widestr('OK'))
    if ischar(casenm)
        % ECLIPSE
        processRestartDiagnostics(casenm, 'outputdir', precompDir);
    elseif isstruct(casenm)
        % If processing MRST. NB casenm = problem for MRST input
        processStatesDiagnostics(casenm, 'outputdir', precompDir);
    end
end
end
