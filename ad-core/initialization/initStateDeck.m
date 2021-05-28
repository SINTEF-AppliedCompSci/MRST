function state = initStateDeck(model, deck)
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

    assert (isstruct(deck) && ...
        all(isfield(deck, {'RUNSPEC', 'SOLUTION'})),   ...
        ['Parameter ''deck'' does not appear to be a ', ...
        'valid ECLIPSE deck']);
    if ~isfield(deck.RUNSPEC, 'SI')
        % We use MRST's gravity to initialize. Must be SI units.
        error('Non-SI units detected. Call ''convertDeckUnits(deck)'' first.')
    end
    if isfield(deck.SOLUTION, 'EQUIL')
        % Equilibrium regions, initialize via solver
        regions = getInitializationRegionsDeck(model, deck);
        [state, p] = initStateBlackOilAD(model, regions);
        if isfield(deck.PROPS, 'SWATINIT')
            model = model.validateModel();
            swat = deck.PROPS.SWATINIT(model.G.cells.indexMap);
            state = model.setProp(state, 'sw', swat);
            so = 1 - swat;
            if model.gas
                so = so - model.getProp(state, 'sg');
            end
            if any(so < 0)
                warning('SWAT resulted in %d negative saturations', sum(so < 0));
            end
            so = max(so, 0);
            state = model.setProp(state, 'so', so);
            pc = model.getProp(state, 'capillarypressure');
            ix = model.getPhaseIndices();
            wix = ix(1);
            oix = ix(2);
            % Note sign
            pcow = -pc{wix};
            if isempty(pcow)
                warning('Oil-Water capillary pressure must be present for scaling of pc-curve with SWATINIT.')
            end
            p(:, oix) = state.pressure; % May be wrong in regions where oil is immobile
            dp = p(:, oix) - p(:, wix);
            bad = dp <= 0;
            scale = dp./pcow;
            scale(bad) = 1;
            state.pcowScale = scale;
        end
    elseif any(isfield(deck.SOLUTION, {'PRESSURE', 'PRVD'})) && ...
            any(isfield(deck.SOLUTION, {'SGAS', 'SOIL', 'SWAT'}))
        % We got the values directly and just need to assign them.
        state = directAssignment(model, deck);
    else
        error(msgid('Scheme:Unknown'), ...
            ['Initialisation scheme specified in ', ...
            '''deck'' is not supported.']);
    end
    if ~isempty(model.AquiferModel)
        state = model.AquiferModel.initStateAquifer(state);
    end
end

function state = directAssignment(model, deck)
    G = model.G;
    phases = model.getPhaseNames();
    nph = numel(phases);
    state = initResSol(model.G, 0);
    state.s = zeros(G.cells.num, nph);
    imap = G.cells.indexMap;
    pvtreg = getPVTNUM(model, deck);
    if isfield(deck.SOLUTION, 'PRVD')
        % Pressure vs depth table
        prvd = deck.SOLUTION.PRVD;
        for reg = 1:numel(prvd)
            F = griddedInterpolant(prvd(:, 1), prvd(:, 2), 'linear', 'nearest');
            act = pvtreg == reg;
            state.pressure(act) = F(G.cells.centroids(act, 3));
        end
    else
        % Pressure specified per cell,
        state.pressure = reshape(deck.SOLUTION.PRESSURE(imap), [], 1);
    end
    % Get saturations, if present
    missing = [];
    for i = 1:nph
        switch lower(phases(i))
            case 'w'
                fld = 'SWAT';
            case 'o'
                fld = 'SOIL';
            case 'g'
                fld = 'SGAS';
        end
        present = isfield(deck.SOLUTION, fld);
        if present
            state.s(:, i) = deck.SOLUTION.(fld)(imap);
        else
            if isempty(missing)
                missing = i;
            else
                error('Too few phase saturations specified');
            end
        end
    end
    if ~isempty(missing)
        state.s(:, missing) = 1 - sum(state.s, 2);
    end
    % Treat RS
    if isprop(model, 'disgas') && model.disgas
        if isfield(deck.SOLUTION, 'RS')
            rs = reshape(deck.SOLUTION.RS(imap), [], 1);
        elseif isfield(deck.SOLUTION, 'PBUB')
            pbub = reshape(deck.SOLUTION.PBUB(imap), [], 1);
            rs = nan(G.cells.num, 1);
            for reg = 1:numel(deck.PROPS.PVTO)
                irs = interp_rs(deck.PROPS.PVTO{reg});
                cells = pvtreg == reg;
                rs(cells) = irs(min(pbub(cells), state.pressure(cells)));
            end
            assert (all(isfinite(rs)), ...
                'Some cells not covered by Solution GOR initialisation.');
        else
            rs = zeros(G.cells.num, 1);
        end
        state.rs = rs;
    end
    % Treat RV
    if isprop(model, 'vapoil') && model.vapoil
        if isfield(deck.SOLUTION, 'RV')
            rv = reshape(deck.SOLUTION.RV(imap), [], 1);
        elseif isfield(deck.SOLUTION, 'PDEW')
            pdew = reshape(deck.SOLUTION.PDEW(imap), [], 1);
            rv = nan(G.cells.num, 1);
            pG = state.pressure;
            if isfield(model.fluid, 'pcOG')
                sG = model.getProp(state, 'sG');
                pG = pG + model.fluid.pcOG(sG);
            end
            
            for reg = 1:numel(deck.PROPS.PVTG)
                irv = interp_rv(deck.PROPS.PVTG{reg});
                cells = pvtreg == reg;
                rv(cells) = irv(min(pdew(cells), pG(cells)));
            end
            assert (all(isfinite(rv)), ...
                'Some cells not covered by Solution OGR initialisation.');
        else
            rv = zeros(G.cells.num, 1);
        end
        state.rv = rv;
    end
    % TODO: Is there a valid direct assigment for these? We should add that.
    if isprop(model, 'polymer')
        state.cp = zeros(G.cells.num, 1);
    end
    if isprop(model, 'surfactant')
        state.cs = zeros(G.cells.num, 1);
    end
end

function reg = getPVTNUM(model, deck)
    if isfield(deck.REGIONS, 'PVTNUM')
        reg = deck.REGIONS.PVTNUM(model.G.cells.indexMap);
    else
        reg = ones(model.G.cells.num, 1);
    end
end

function irs = interp_rs(PVTO)
   irs = @(p) interp1(PVTO.data(PVTO.pos(1:end-1), 1), ...
                      PVTO.key                       , ...
                      p, 'linear', 'extrap');
end
