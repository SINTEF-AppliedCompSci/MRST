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
%   trapAnalysis, showTrappingStructure

%{
#COPYRIGHT#
%}
    ff = fopen(filename);
    meta = struct();
    while 1
        line = fgetl(ff);
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
