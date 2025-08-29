function flatjsonviewer = flattenJsonStruct(jsonstruct, varargin)
%
% A json struct is hierarchical by design. For visualization, it is however convenient to have a flattened version with
% all the entries given in the structure at the same level. This function returns an object of the class
% :battmo:`FlatJsonViewer` which offers visualization capabilities with sorting and filtering.
    
    opt = struct('doprint', true);
    opt = merge_options(opt, varargin{:});
    
    flatjson = flattenJsonStruct_({}, jsonstruct, []);
    flatjson = reshape(flatjson, 2, [])';

    flatjsonviewer = FlatJsonViewer(flatjson);

    if opt.doprint
        flatjsonviewer.print();
    end

end

function flatjson = flattenJsonStruct_(flatjson, jsonstruct, prefix)

    dostruct = false;
    
    if isobject(jsonstruct)
        fds = properties(jsonstruct);
        dostruct = true;
    elseif isstruct(jsonstruct)    
        fds = fieldnames(jsonstruct);
        dostruct = true;
    else
        flatjson{end + 1} = prefix;
        flatjson{end + 1} = jsonstruct;
    end

    if dostruct
        for ifd = 1 : numel(fds)
            fd = fds{ifd};
            if isempty(prefix)
                subprefix = fd;
            else
                subprefix = sprintf('%s.%s', prefix, fd);
            end
            subjsonstruct = jsonstruct.(fd);
            if iscell(subjsonstruct) && numel(subjsonstruct) > 1
                subjsonstruct = {subjsonstruct};
            end
            flatjson = flattenJsonStruct_(flatjson, subjsonstruct, subprefix);
        end
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
