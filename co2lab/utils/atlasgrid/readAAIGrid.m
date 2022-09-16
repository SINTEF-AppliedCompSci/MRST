function [meta, data] = readAAIGrid(filename)
% Read AIIGrid from file.
%
% SYNOPSIS:
%       [meta, data] = readAAIGrid('~/testfile')
%
% PARAMETERS:
%   inp     - A valid filename for fopen.
%
% RETURNS:
%   meta    - metadata defining global coordinate system, cell size, etc.
%   
%   data    - A matrix containing the data from the file.
%
%
% SEE ALSO:
%   `trapAnalysis`, `showTrappingStructure`

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    [ff, msg] = fopen(filename, 'rt');
    if ff < 0
        error('File:Open', 'Failed to open ''%s'': %s', filename, msg);
    end

    cleanup = onCleanup(@() fclose(ff));

    meta = struct();
    while 1
        line = fgetl(ff);
        
        % replace any commas present with decimal points (which are present
        % in some of the Norwegian Sea formation data)
        line = strrep(line,',','.');
        
        subs = regexp(line, ' *','split');
        if numel(subs{1}) == 0 || numel(subs) > 2
            break
        end
        meta.(subs{1}) = sscanf(subs{2}, '%e');
    end
    
    data = zeros([meta.nrows, meta.ncols]);
    i = 1;
    while i < meta.nrows;
        line = fgetl(ff);
        
        % replace any commas present with decimal points (which are present
        % in some of the Norwegian Sea formation data)
        line = strrep(line,',','.');
        
        vals = sscanf(line, '%e');
        data(i, 1:numel(vals)) = vals;
        i = i + 1;
    end
    data(data == meta.NODATA_value) = nan;
    data(isinf(data))               = nan;
    
    % Mask zero values to not a number
    data(data == 0) = nan;
    
    meta.dims = [meta.ncols, meta.nrows];
    
    % Reverse to get deck files with the same result as johansen/utsira
    % datasets.
    data(:,:) = data(end:-1:1,:);
    
    % Switch from C style indexing to FORTRAN/MATLAB-style
    data = data'; 

end
