function viewStruct(mstruct, input)

    if nargin < 2
        input = UnAssigned;
    end

    ca = getStructField(input, 'clean_up_large_arrays', false);

    if ca
        n  = getStructField(input, 'array_size', 10);
        mstruct = clean_up_large_arrays(mstruct, n);
    end

    rfds = getStructField(input, 'remove_fields', {}); % cell array of strings

    for irfd = 1 : numel(rfds)
        fd = rfds{irfd};
        mstruct = setStructField(mstruct, fd, 'removed for printing', 'handleMisMatch', 'quiet');
    end
    
    display(jsonencode(mstruct, 'PrettyPrint', true));

end

function cleaned = clean_up_large_arrays(mstruct, n)

    if isstruct(mstruct)
        fields = fieldnames(mstruct);
        for i = 1:length(fields)
            field = fields{i};
            mstruct.(field) = clean_up_large_arrays(mstruct.(field), n);
        end
        cleaned = mstruct;

    elseif iscell(mstruct)
        for i = 1:length(mstruct)
            mstruct{i} = clean_up_large_arrays(mstruct{i}, n);
        end
        cleaned = mstruct;

    elseif isnumeric(mstruct) || islogical(mstruct) || isstring(mstruct) || ischar(mstruct)
        if numel(mstruct) > n
            cleaned = sprintf('Array of size %s', mat2str(size(mstruct)));
        else
            cleaned = mstruct;
        end
    else
        cleaned = mstruct;
    end

end


%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
