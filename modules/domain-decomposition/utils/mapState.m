function state = mapState(state, substate, mappings, varargin)
% Map substate to a full state

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

    opt = struct('includeOverlap' , false, ...
                 'includeExternal', false, ...
                 'mapCellFields'  , true , ...
                 'mapFaceFields'  , false, ...
                 'mapWellSol'     , true );
    opt = merge_options(opt, varargin{:});
    if opt.mapCellFields
        % Map cell fields
        glob_cells = mappings.cells.internal;
        if opt.includeOverlap
            glob_cells = glob_cells | mappings.cells.overlap;
        end
        if opt.includeExternal
            glob_cells = glob_cells | mappings.cells.external;
        end
        loc_cells = glob_cells(mappings.cells.keep);
        fields = union(mappings.cells.fields, {'dpRel', 'dpAbs'});
        state = mapFields(state, substate, fields, glob_cells, loc_cells);
    end
    if opt.mapFaceFields
        % Map face fields
        glob_faces = mappings.faces.internal;
        if opt.includeOverlap
            glob_faces = glob_faces | mappings.faces.overlap;
        end
        if opt.includeExternal
            glob_faces = glob_faces | mappings.faces.external;
        end
        loc_faces = glob_faces(mappings.faces.keep);
        fields = mappings.faces.fields;
        state = mapFields(state, substate, fields, glob_faces, loc_faces);
    end
    if opt.mapWellSol && ~isempty(substate.wellSol)
        % Map well solution
        wfields = fieldnames(substate.wellSol);
        state.wellSol = mapWellSol(state.wellSol, substate.wellSol, wfields, mappings.wells.keep);
    end
end

%-------------------------------------------------------------------------%
function state = mapFields(state, substate, fields, ix, subix)
    for j = 1:numel(fields)
        f = fields{j};
        if ~isfield(substate, f)
            continue
        end
        if ~isfield(state, f)
            state.(f) = zeros(numel(ix), size(substate.(f), 2));
        end
        if size(state.(f), 2) < size(substate.(f), 2)
            v = state.(f);
            state.(f) = zeros(size(state.(f),1), size(substate.(f),2));
            state.(f)(:,1) = v;
        end
        state.(f)(ix,:) = substate.(f)(subix,:);
    end
end

%-------------------------------------------------------------------------%
function wellSol = mapWellSol(wellSol, subWellSol, fields, ix)
    for j = 1:numel(fields)
        f = fields{j};
        if ~isfield(wellSol, f)
            [wellSol.(f)] = deal([]);
        end
        [wellSol(ix).(f)] = subWellSol.(f);
    end
end
