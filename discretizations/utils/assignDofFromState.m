function state = assignDofFromState(disc, state, names)
    % Assign dofs from state (typically initial state). All dofs
    % for dofNo > 0 are set to zero.
    % Set degree to disc.degree in all cells
    if ~isfield(state, 'degree')
        state.degree = repmat(disc.degree, disc.G.cells.num, 1);
    end
    % Create vector dofPos for position of dofs in state.sdof
    if 0
        state = disc.updateDofPos(state);
        state.nDof  = disc.getnDof(state);
    else
        state.nDof = repmat(disc.basis.nDof, disc.G.cells.num, 1);
        state.nDof(all(state.degree == 0,2)) = 1;
        state = disc.updateDofPos(state);
    end
    ix          = disc.getDofIx(state, 1);
    % Loop trough possible fields and initialise constant dof
    if nargin < 3
        names = fieldnames(state);
    end
    except = exceptions();
    for fNo = 1:numel(names)
        f = names{fNo};
        if isfield(state, f) && ~any(strcmpi(f, except))
            if size(state.(f),1) == disc.G.cells.num
                dof       = zeros(sum(state.nDof), size(state.(f),2));
                dof(ix,:) = state.(f);
                state.([f, 'dof']) = dof;
            else
                v = repmat(state.(f), sum(state.nDof),1);
                state.([f, 'dof']) = v;
            end
        end
    end

end

function flds = exceptions()
    flds = {'nDof', 'degree', 'wellSol', 'flux', 'sMax', 'dofPos', 'eos', 'flag', 'FlowProps'};
end