function state = assignDofFromState(disc, state, names)
    % Assign dofs from state (typically initial state). All dofs
    % for dofNo > 0 are set to zero.
    % Set degree to disc.degree in all cells

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
