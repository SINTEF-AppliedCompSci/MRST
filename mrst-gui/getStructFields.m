function flds = getStructFields(G, s, name)
%Dump possible plotting fields of a struct into human readable format.
%
% SYNOPSIS:
%   flds = getStructFields(G, rock, 'dataset')
%
% DESCRIPTION:
%   getStructFields finds all plottable struct fields (typically
%   G.cells.num long vectors) and returns them as a set of possible strings
%   to be passed onto readStructField.
%
% REQUIRED PARAMETERS:
%   G     - Target grid for plotting.
%
%   s     - Struct with fields suitable for plotting, or simply a numerical
%           array.
%
% RETURNS:
%   A list of names in the form:
%          'rock.poro', 'rock.perm:1', 'rock.perm:2' etc.
%
% SEE ALSO:
%   readStructFields, datasetSelector

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

    flds = {};
    if isnumeric(s)
        if size(s, 1) == G.cells.num
            N = min(size(s, 2), 1000);
            if N > 1
                nn = arrayfun(@(x) [name, ':', num2str(x)], (1:N).', ...
                    'UniformOutput', false);

                flds = [flds; nn];
            end

            if N == 3 || N == 1
                % special case true color or single dataset
                flds = [flds; name];
            end
        end
        return
    end
    if isstruct(s)
        f = fieldnames(s);
    else
        f = {};
    end

    if ~isempty(s)
        for i = 1:numel(f)
            subs = s.(f{i});
            next = getStructFields(G, subs, f{i});
            if ~isempty(next)
                if ischar(next{1})
                    flds = [flds; strcat([name '.'], next)];
                end
            end
        end
    end
end
