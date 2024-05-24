function wellSols = getWellSolsFromWellboreState(model, states, wbname)
%Extract the wellbore states from a set of states of a simulation with a CompositeModel where one of the models is a WellboreModel
%
% SYNOPSIS:
%   wellSols = getWellSolsFromWellboreState(model, states)
%   wellSols = getWellSolsFromWellboreState(model, states, wbname)
%
% DESCRIPTION:
%   This function extracts properties from state.Wellbore for each state
%   in states and converts it to a the standard wellSol format of MRST. The
%   resulting wellSols can be used in e.g., plotWellSols.
%
% REQUIRED PARAMETERS:
%   model  - Instance of CompositeModel where one of the submodels
%            is a WellboreModel
%
%   states - States from a simulation with the model
%
% OPTIONAL PARAMETERS:
%   wbname - Name of the wellbore model as defined in
%            model.submodels.(wbname). Default: Wellbore
%
% RETURNS:
%   wellSol - Well solution structure array
%
% NOTE:
%   The WellboreModel is currently only compatible with GeothermalModel

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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
    
    % Set wellbore name if not given
    if nargin < 3, wbname = 'Wellbore'; end
    
    % Field to output in each wellSol
    names = {'bhp', 'bht', 'effect', 'qWs', 'qOs', 'qGs'};
    
    n = numel(states);
    if isa(states, 'ResultHandler'), n = states.numelData(); end

    wellSols   = cell(n,1);
    wellNames  = {model.submodels.(wbname).wells.name};
    groupNames = {};
    if model.submodels.(wbname).numGroups > 0
        groupNames = {model.submodels.(wbname).groups.name};
    end
    
    nw = model.submodels.(wbname).numWells();
    ng = model.submodels.(wbname).numGroups();
    v0 = cell(1,numel(names)+1);
    
    for i = 1:n
        state = states{i};
        state = state.(wbname);
        ws = repmat(cell2struct(v0, ['name', names], 2), 1, nw + ng);
        [ws.name] = deal(wellNames{:}, groupNames{:});
        for name = names
            if ~isfield(state, name{1}), [ws.(name{1})] = deal([]); continue; end
            v = model.submodels.(wbname).getProp(state, name{1});
            if strcmpi(name, 'bht'), name = {'T'}; end %#ok
            v = num2cell(v);
            [ws.(name{1})] = deal(v{:});
        end
        wellSols{i} = ws;
    end
    
end