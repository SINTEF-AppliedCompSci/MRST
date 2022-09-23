function PROPS = generatePROPS(model, varargin)
% Generate PROPS from model. Very limited function.

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

    opt = struct('pMax', 500*barsa,...
                 'pMin', 1*atm, ...
                 'ns', 21,...
                 'np', 21,...
                 'props', struct(), ...
                 'pRef', 50*barsa, ...
                 'writeExtra', false);
    opt = merge_options(opt, varargin{:});
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
    PROPS = opt.props;
    if ~isfield(PROPS, 'DENSITY')
        PROPS.DENSITY = density;
    end

    if model.water && ~(isfield(PROPS, 'PVTW') || isfield(PROPS, 'PVTW_EXTENDED'))
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
            if ~isfield(PROPS, 'PVTG')
                warning('Writing PVTG from fluid object is not supported.')
            end
        elseif ~is_comp && ~isfield(PROPS, 'PVDG')
            % PVT table
            bg = f.bG(p);
            mug = f.muG(p);
            PROPS.PVDG = {[p, 1./bg, mug]};
        end
        % Saturation table
        if ~isfield(PROPS, 'SGOF')
            krg = f.krG(s);
            if isfield(f, 'krOG')
                krog = f.krOG(1-s);
            else
                krog = f.krO(1-s);
            end
            if isfield(f, 'pcOG')
                pc = f.pcOG(s);
            else
                pc = zeros(opt.ns, 1);
            end
            PROPS.SGOF = {[s, krg, krog, -pc]};
        end
    end

    if model.water && model.oil
        % PVT table
        if isprop(model, 'disgas') && model.disgas
            if ~isfield(PROPS, 'PVTO')
                warning('Writing PVTO from fluid object is not supported.')
            end
        elseif ~is_comp && ~isfield(PROPS, 'PVDO')
            bo = f.bO(p);
            muo = f.muO(p);
            PROPS.PVDO = {[p, 1./bo, muo]};
        end
        if ~isfield(PROPS, 'SWOF')
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
    end
    % TODO: pvMultR...
    if isfield(f, 'pvMultR')
        warning('pvMultR not supported... Ignored.');
    end
    if ~isfield(PROPS, 'ROCK')
        PROPS.ROCK = [1*atm, 0, NaN, NaN, NaN, NaN];
    end
end
