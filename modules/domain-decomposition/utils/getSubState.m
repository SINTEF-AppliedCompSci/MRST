function [substate, mappings] = getSubState(state, mappings, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('mapCellFields'    , true , ...
                 'mapFaceFields'    , true , ...
                 'mapWellSol'       , true , ...
                 'mapStateFunctions', true, ...
                 'useMatrix'        , false);
    opt = merge_options(opt, varargin{:});
    % Make copy
    substate = state;
    % Map cell fields
    mappings.cells.fields = {};
    if opt.mapCellFields
        [substate, cellFields] = mapFields(state, substate, mappings.cells, opt);
        mappings.cells.fields  = cellFields;
    end
    % Map face fields
    mappings.faces.fields = {};
    if opt.mapFaceFields
        [substate, faceFields] = mapFields(state, substate, mappings.faces, opt);
        mappings.faces.fields  = faceFields;
    end
    % Map well solution
    if opt.mapWellSol && isfield(state, 'wellSol') && isfield(mappings, 'wells')
        substate.wellSol = state.wellSol(mappings.wells.keep);
        if isfield(state, 'FacilityState')
            wmap = mappings.wells;
            wmap.keep = wmap.keep([state.wellSol.status]);
            substate.FacilityState ...
                = mapFields(state.FacilityState, substate.FacilityState, wmap, opt);
        end
    end
    % Map state functions
    if opt.mapStateFunctions
        if isfield(state, 'FlowProps')
            substate.FlowProps = mapFields(state.FlowProps, substate.FlowProps, mappings.cells, opt);
        end
        if isfield(state, 'PVTProps')
            substate.PVTProps = mapFields(state.PVTProps, substate.PVTProps, mappings.cells, opt);
        end
        if isfield(state, 'FluxProps')
            substate.FluxProps = mapFields(state.FluxProps, substate.FluxProps, mappings.cells, opt);
            substate.FluxProps = mapFields(state.FluxProps, substate.FluxProps, mappings.faces, opt);
        end
    end
    % Map pressure reduction factors
    if isfield(state, 'reductionFactorProps')
        substate.reductionFactorProps = mapFields(state.reductionFactorProps, substate.reductionFactorProps, mappings.cells, opt);
    end

end

%-------------------------------------------------------------------------%
function [substate, fields] = mapFields(state, substate, map, opt)
    n = numel(map.keep);
    if isa(state, 'HandleStruct')
        fields = fieldnames(state.struct);
    else
        fields = fieldnames(state);
    end
    keep   = false(numel(fields), 1);
    if opt.useMatrix
        M = sparse(1:nnz(map.keep), find(map.keep), 1, nnz(map.keep), numel(map.keep));
        getSubset = @(v) M*v;
    else
        getSubset = @(v) v(map.keep,:);
    end
    get = @(v) (isnumeric(value(v)) || islogical(v)) && size(value(v), 1) == n;
    for i = 1:numel(fields)
        f = fields{i};
        v = state.(f); 
        if ~iscell(v)
            if get(v)
                keep(i)      = true;
                substate.(f) = getSubset(v);
            end
        else
            for j = 1:numel(v)
                if get(v{j})
                    keep(i)         = true;
                    substate.(f){j} = getSubset(v{j});
                end
            end
        end
    end
    fields = fields(keep);
end
