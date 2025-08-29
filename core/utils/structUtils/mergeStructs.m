function mstruct = mergeStructs(mstructs, varargin)
%
%  We call a json structure, abbreviated mstruct, a MATLAB structure produced by the jsondecode command or :battmo:`parseBattmoJson`
%    
%  The input mstructs is a list of mstruct
%    
%  The command mergeStructs merges recursively all the mstruct contained in the list mstructs
%    
%  If two mstruct assign the same field, then the first one is applied and a warning message is sent. use warn=false to switch off this message

    opt = struct('warn', true);
    opt.tree = {};
    opt = merge_options(opt, varargin{:});

    % needed to produce meaningfull error message
    tree = opt.tree;
    
    if numel(mstructs) == 1
        mstruct = mstructs{1};
        return;
    end
    
    mstruct1 = mstructs{1};
    mstruct2 = mstructs{2};
    if numel(mstructs) >= 2
        mstructrests = mstructs(3 : end);
    else
        mstructrests = {};
    end

    mstruct = mstruct1;

    if isstruct(mstruct) && numel(mstruct) > 1
        % We look at the case where the mstruct is a struct-array

        nelts = numel(mstruct);

        if nelts ~= numel(mstruct2)

            % we overwrite the elements, which means do nothing
            return
            
        else
            % we merge each of the elements

            for ielt = 1 : nelts
                subtree = tree;
                subtree{end + 1} = sprintf('[%d]', ielt);
                mstruct(ielt) = mergeStructs({mstruct(ielt), mstruct2(ielt)}, 'warn', opt.warn, 'tree', subtree);
            end

            return
            
        end

    end
    
    fds1 = fieldnames(mstruct1);
    fds2 = fieldnames(mstruct2);
    
    for ifd2 = 1 : numel(fds2)
        
        fd2 = fds2{ifd2};
        
        if ~ismember(fd2, fds1)
            % ok, add the substructure
            mstruct.(fd2) = mstruct2.(fd2);
        else
            if isstruct(mstruct.(fd2)) && isstruct(mstruct2.(fd2))
                % we have to check the substructure
                subtree = tree;
                subtree{end + 1} = fd2;
                submstruct = mergeStructs({mstruct.(fd2), mstruct2.(fd2)}, 'warn', opt.warn, 'tree', subtree);
                mstruct.(fd2) = submstruct;
            elseif ~isstruct(mstruct.(fd2)) && ~isstruct(mstruct2.(fd2)) && isequal(mstruct.(fd2), mstruct2.(fd2))
                % ok. Both are given but same values
            elseif opt.warn
                varname = horzcat(tree, fd2);
                fprintf('mergeStructs: Parameter %s is assigned twice with different values. Value from first mstruct is used.\n', strjoin(varname, '.'));
            end
        end
    end

    if ~isempty(mstructrests)
        mstruct = mergeStructs({mstruct, mstructrests{:}}, 'warn', opt.warn, 'tree', {});
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
