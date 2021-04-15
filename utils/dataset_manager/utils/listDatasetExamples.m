function varargout = listDatasetExamples(name)
% List all MRST examples using a specific dataset
%
% SYNOPSIS:
%   listDatasetExamples('mydataset');
%   ex = listDatasetExamples('mydataset');
%
% DESCRIPTION:
%   This function has two calling syntaxes. If called without output
%   arguments, it will print a list of examples using a given example into
%   the command line, with links to the corresponding files.
%
%   If called with a single output argument, it will instead output the
%   list of examples as a variable, while not printing anything to the
%   command window.
%
% REQUIRED PARAMETERS:
%   name   - Either a string containing a valid dataset name, or a info
%            struct for a given dataset.
%
% RETURNS:
%   examples - (OPTIONAL) A array of structs, each representing an example
%              where the dataset is used. Contains fields name, path and
%              module.
%
% NOTE:
%    This function relies on the `dataset_*.m` functions being up to date
%    with the examples included in MRST. Although we strive to do so, it
%    may happen that examples are not listed.
%
% SEE ALSO:
%   `mrstDatasetGUI`

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


    if isstruct(name)
        info = name;
    elseif ischar(name)
        info = getDatasetInfo(name);
    else
        error('Invalid input for dataset name');
    end
    
    nex = numel(info.examples);
    for i = 1:nex
        ex = infostruct(info.examples{i});
        if i == 1
            examples = repmat(ex, nex, 1);
        else
            examples(i) = ex;
        end
    end
    if nex==0, examples = {}; end
    
    if nargout > 0
        varargout{1} = examples;
    else
        printExampleList(info, examples);
    end
end

function s = infostruct(name)
    v = regexp(name, ':', 'split');
    if numel(v) == 1
        % Example is a part of core
        module = '';
        fn = v{1};
    else
        % Example is a part of a module.
        module = v{1};
        fn = v{2};
    end
    
    if ~isempty(module)
        prevMod = mrstModule();
        % Ensure module is in path, and was the most recently added to
        % avoid confusion with examples that have the same name
        mrstModule('add', module);
        pth = which(fn);
        mrstModule('reset', prevMod{:});
    else
        pth = which(fn);
    end
    s = struct('name', fn, 'module', module, 'path', pth);
end

function printExampleList(info, examples)
    N = max(arrayfun(@(x) numel(x.name), examples));

    fprintf('\n*** List of examples using ''%s'' dataset: ***\n', info.name)
    
    for i = 1:numel(examples)
        ex = examples(i);
        fprintf(['%-', num2str(N), 's (%s): <a href="matlab: mrstModule add %s, edit ''%s''">%s</a>\n'], ...
                ex.name, ex.module, ex.module, ex.path, ex.path);

    end
end
