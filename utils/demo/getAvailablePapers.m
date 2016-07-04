function papers = getAvailablePapers()
%
%
% SYNOPSIS:
%   
%
% DESCRIPTION:
%   
%
% REQUIRED PARAMETERS:
%
% RETURNS:
%   
% SEE ALSO:
%   

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
    pth = mfilename('fullpath');
    pth = pth(1:(end-numel(mfilename()) - 1));
    sets = dir(fullfile(pth, 'papers'));
    
    % Extract functions that start with dataset and is not directories
    sn = 'paper';
    ok = strncmp({sets.name}, sn, numel(sn)) & ~[sets.isdir];

    names = {sets(ok).name}';
    
    nD = numel(names);
    
    % Loop over candidate functions, adding to array as we go (dynamic
    % expansion is ok since the number of datasets is relatively small).
    papers = [];
    for i = 1:nD
        [name, name, ext] = fileparts(names{i}); %#ok
        if ~strcmpi(ext, '.m')
            continue
        end
        
        paper = eval([name, '()']);
        papers = [papers; paper];
    end
end
