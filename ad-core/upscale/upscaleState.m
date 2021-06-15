function state = upscaleState(coarsemodel, model, state)
%Create a upscaled state by simple processing of values
%
% SYNOPSIS:
%   state_coarse = upscaleState(coarsemodel, model, state_fine)
%
% DESCRIPTION:
%   Convert a state for a fine model into a realization of the same state
%   for a coarse model.
%
% REQUIRED PARAMETERS:
%   coarsemodel - A coarse model derived from the fine model.
%
%   model       - The fine model. Subclass of `ReservoirModel`.
%
%   state       - State to be converted. Should correspond to `model`.
%
%
% RETURNS:
%   state       - Coarse state suitable for the coarsemodel.
%
% SEE ALSO:
%   `upscaleSchedule`, `upscaleModelTPFA`

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
    p = coarsemodel.G.partition;
    CG = coarsemodel.G;
    
    pvc = coarsemodel.operators.pv;
    pvf = model.operators.pv;
    
    counts = accumarray(p, 1);
    state_f = state;
    
    nph = size(state.s, 2);
    pvs = bsxfun(@times, state.s, pvf);
    
    % Calculate saturations based on new and old pore volume
    s = zeros(CG.cells.num, nph);
    for i = 1:nph
        s(:, i) = accumarray(p, pvs(:, i))./pvc;
    end
    state.s = s;
    if isfield(state, 'components')
        rhoL = model.PropertyModel.computeMolarDensity(model.EOSModel, state.pressure, state_f.x, state_f.Z_L, state.T, true);
        rhoV = model.PropertyModel.computeMolarDensity(model.EOSModel, state.pressure, state_f.y, state_f.Z_V, state.T, false);
        % Get saturations for hydrocarbon-like phases
        [sL, sV] = model.getProps(state_f, 'so', 'sg');
        % Calculate total molar density
        N = rhoL.*state_f.x.*sL + rhoV.*state_f.y.*sV;
        % Fine number of moles
        N_f = bsxfun(@times, N, pvf);
        ncomp = size(state.components, 2);
        % Compute coarse moles
        N_c = zeros(CG.cells.num, ncomp);
        % Coarse molar density
        for i = 1:ncomp
            N_c(:, i) = accumarray(p, N_f(:, i))./pvc;
        end
        state.components = bsxfun(@rdivide, N_c, sum(N_c, 2));
        flds = {'L', 'x', 'y', 'K', 'K', 'Z_L', 'Z_V', 'mixing', 'flag', 'eos'};
        for i = 1:numel(flds)
            f = flds{i};
            if isfield(state, f)
                state = rmfield(state, f);
            end
        end
    end
    
    if isprop(model, 'disgas') && model.disgas && isfield(state, 'rs')
        sO = model.getProp(state_f, 'sO');
        sO_c = coarsemodel.getProp(state, 'sO');
        state.rs = accumarray(p, sO.*state.rs.*pvf)./(sO_c.*pvc);
        state.rs(~isfinite(state.rs)) = 0;
    end
    if isprop(model, 'vapoil') && model.vapoil && isfield(state, 'rv')
        sG = model.getProp(state_f, 'sG');
        sG_c = coarsemodel.getProp(state, 'sG');
        state.rv = accumarray(p, sG.*state.rv.*pvf)./(sG_c.*pvc);
        state.rv(~isfinite(state.rv)) = 0;
    end
    % Average the pressure (not entirely correct for compressible systems,
    % but we won't start evaluating properties in here).
    %state.pressure = accumarray(p, state.pressure)./counts;
    state.pressure = accumarray(p, pvf.*state.pressure)./pvc;
    if isfield(state, 'T')
        state.T = accumarray(p, state.T)./counts;
    end
    if isfield(state, 'flux')
        cfsign = fineToCoarseSign(CG);
        cfacesno = rldecode(1:CG.faces.num, diff(CG.faces.connPos), 2) .';
        newflux = zeros(CG.faces.num, size(state.flux, 2));
        for i = 1:size(state.flux, 2)
            newflux(:, i)   = accumarray(cfacesno, state.flux(CG.faces.fconn, i) .* cfsign);
        end
        state.flux = newflux;
    end
    if isfield(state, 'FlowProps')
        state = rmfield(state, 'FlowProps');
    end
    if isfield(state, 'components')
        % Flash compositional variables
        state = model.computeFlash(state, inf);
    end
end
