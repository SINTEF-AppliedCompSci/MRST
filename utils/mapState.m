function state = mapState(state, substate, mappings, varargin)
    % Map substate to a full state
    opt = struct('includeOverlap', false, ...
                 'mapCellFields' , true , ...
                 'mapFaceFields' , false, ...
                 'mapWellSol'    , true );
    opt = merge_options(opt, varargin{:});
    if opt.mapCellFields
        % Map cell fields
        glob_cells = mappings.cells.internal;
        if opt.includeOverlap
            glob_cells = glob_cells | mappings.cells.overlap;
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