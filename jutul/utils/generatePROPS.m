function PROPS = generatePROPS(model, varargin)
    opt = struct('pMax', 500*barsa,...
                 'pMin', 1*atm, ...
                 'ns', 21,...
                 'np', 21,...
                 'pRef', 50*barsa, ...
                 'writeExtra', false);
    opt = merge_options(opt, varargin{:});
    if isprop(model, 'disgas') && (model.disgas || model.vapoil)
        error('Creating props only supported for immiscibile (deadoil) case.');
    end
    is_comp = isa(model, 'ThreePhaseCompositionalModel');
    pR = opt.pRef;
    density = [1, 1, 1];
    f = model.fluid;
    p = linspace(opt.pMin, opt.pMax, opt.np)';
    s = linspace(0, 1, opt.ns)';
    if model.water
        density(2) = f.rhoWS;
    end
    if model.oil
        density(1) = f.rhoOS;
    end
    if model.gas
        density(3) = f.rhoGS;
    end
    PROPS = struct();
    PROPS.DENSITY = density;

    if model.water
        % PVT table
        write_pvtw = true;
        if opt.writeExtra
            bw = f.bW(p);
            muw = f.muW(p);
            % If viscosity or b depends on pressure at all, add that table
            if bw(1) ~= bw(end) || muw(1) ~= muw(end)
                PROPS.PVTW_EXTENDED = {[p, 1./bw, muw]};
                write_pvtw = false;
            end
        end
        if write_pvtw
            bWr = f.bW(pR);
            mur = f.muW(pR);
            PROPS.PVTW = [pR, bWr, 0, mur, 0];
        end
    end

    if model.gas && model.oil
        if isprop(model, 'vapoil') && model.vapoil
            warning('VAPOIL not supported for PVT table.')
        elseif ~is_comp
            % PVT table
            bg = f.bG(p);
            mug = f.muG(p);
            PROPS.PVDG = {[p, 1./bg, mug]};
        end
        % Saturation table
        krg = f.krG(s);
        if isfield(f, 'krOG')
            krog = f.krO(1-s);
        else
            krog = f.krOG(1-s);
        end
        if isfield(f, 'pcOG')
            pc = f.pcOG(s);
        else
            pc = zeros(opt.ns, 1);
        end
        PROPS.SGOF = {[s, krg, krog, -pc]};
    end

    if model.water && model.oil
        % PVT table
        if isprop(model, 'disgas') && model.disgas
            warning('DISGAS not supported for PVT table.')
        elseif ~is_comp
            bo = f.bO(p);
            muo = f.muO(p);
            PROPS.PVDO = {[p, 1./bo, muo]};
        end
        % Saturation table
        krw = f.krW(s);
        if isfield(f, 'krO')
            krow = f.krO(1-s);
        else
            krow = f.krOW(1-s);
        end
        if isfield(f, 'pcOW')
            pc = f.pcOW(s);
        else
            pc = zeros(opt.ns, 1);
        end
        PROPS.SWOF = {[s, krw, krow, pc]};
    end
    % TODO: pvMultR...
    if isfield(f, 'pvMultR')
        warning('pvMultR not supported... Ignored.');
    end
    PROPS.ROCK = [1*atm, 0, NaN, NaN, NaN, NaN];
end
