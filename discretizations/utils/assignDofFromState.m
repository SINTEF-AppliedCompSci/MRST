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
    flds = getDofFields();
    for f = flds
        if isfield(state, f)
            dof       = zeros(sum(state.nDof), size(state.(f{1}),2));
            dof(ix,:) = state.(f{1});
            state.([f{1}, 'dof']) = dof;
        end
    end

end

function flds = getDofFields()
    flds = {'s', 'rs', 'rv', 'x', 'y', 'components'};
end