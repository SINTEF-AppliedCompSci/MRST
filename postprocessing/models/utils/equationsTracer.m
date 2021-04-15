function [problem, state] = equationsTracer(state0, state, model, dt, drivingForces, varargin)
%Undocumented Utility Function

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

    nt = model.getNumberOfTracers();
    op = model.operators;
    
    [tracers, tracers0] = deal(cell(nt, 1));
    primaryVars = model.tracerNames;
    [tracers{:}] = model.getProps(state, primaryVars{:}, 'wellSol');
    [tracers0{:}] = model.getProps(state0, primaryVars{:});
    
    [tracers{:}] = initVariablesADI(tracers{:});
    %[vT, qT] = getFlux(model, state); 
    vT = sum(state.flux(op.internalConn, :), 2);
    
    [eqs, types] = deal(cell(1, nt));
    names = primaryVars;
    flag = vT > 0;
    
    W = drivingForces.W;
    if ~isempty(W)
        aw = find([W.status]);
        qT = sum(vertcat(state.wellSol(aw).flux), 2);
        wc = vertcat(W(aw).cells);
        p2w = getPerforationToWellMapping(W(aw));
        comp = vertcat(W(aw).tracer);
        comp = comp(p2w, :);
        cmap = sparse(wc, 1:numel(wc), 1, model.G.cells.num, numel(wc));
        inj = qT > 0;
    end

    bc = drivingForces.bc;
    % bc term requires following fields
    % - flux2qextmap : mapping from flux values to flux on external faces
    %                  oriented towards outer domain
    % - flux2qbcmap : mapping from flux values to flux on the source face
    %                  oriented towards outer domain
    % - bc2extmap : mapping from values at source faces to external faces
    % - ext2cellmap : mapping from values at external faces to adjacent cell
    %                 (used to compute influx)
    % - tracers : concentration of the tracers at the source faces.
    
    if ~isempty(bc)
        flux2qextmap = bc.flux2qextmap;
        flux2qbcmap  = bc.flux2qbcmap;
        bc2extmap    = bc.bc2extmap;
        qext2cellmap = bc.qext2cellmap;
        qext = flux2qextmap*state.flux;
        qext(qext < 0) = 0;
        qbc = flux2qbcmap*state.flux;
        qbc(qbc > 0) = 0;
    end

    for i = 1 : nt
        t = tracers{i};
        t0 = tracers0{i};
        vi = op.faceUpstr(flag, t).*vT;
        % Conservation eqn
        eqs{i} = (op.pv./dt).*(t - t0) + op.Div(vi);
        
        if ~isempty(W)
            % Composition source terms
            qi = (inj.*comp(:, i) + ~inj.*t(wc)).*qT;
            eqs{i} = eqs{i} - cmap*qi;
        end
        
        if ~isempty(bc)
            % Composition source terms
            conc = qext2cellmap'*t;
            bcconc = bc.tracers{i};
            q = qext2cellmap*(conc.*qext + bc2extmap*(bcconc.*qbc));
            eqs{i} = eqs{i} + q;
        end
        
        types{i} = 'cell';
    end
    
    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
end

%function [vT, qT] = getFlux(model, state)
%    vT = sum(state.flux(model.operators.internalConn, :), 2);
%    qT = sum(vertcat(state.wellSol.flux), 2);
%end
