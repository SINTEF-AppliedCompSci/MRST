function state = assignDofFromState(disc, state)
    % Assign dofs from state (typically initial state). All dofs
    % for dofNo > 0 are set to zero.

    % Set degree to disc.degree in all cells
    state.degree = repmat(disc.degree, disc.G.cells.num, 1);

    % Create vector dofPos for position of dofs in state.sdof
    state       = disc.updateDofPos(state);
    state.nDof  = disc.getnDof(state);
    ix          = disc.getDofIx(state, 1);
    % Loop trough possible fields and initialise constant dof
%     flds = getDofFields();
    flds   = fieldnames(state);
    except = exceptions();
    for fNo = 1:numel(flds)
        f = flds{fNo};
        if isfield(state, f)                ...
                && ~any(strcmpi(f, except)) ...
                && size(state.(f),1) == disc.G.cells.num
            dof       = zeros(sum(state.nDof), size(state.(f),2));
            dof(ix,:) = state.(f);
            state.([f, 'dof']) = dof;
        end
    end

end

function flds = exceptions()
    flds = {'nDof', 'degree', 'wellSol', 'flux'};
end

% function flds = getDofFields()
%     flds = {'pressure', 's', 'rs', 'rv', 'x', 'y', 'components', 'c'};
% end