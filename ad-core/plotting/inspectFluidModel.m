function varargout = inspectFluidModel(model, varargin)
% Create a interactive plotting panel for a given model that shows
% different fluids properties.
%
% SYNOPSIS:
%   h = inspectFluidModel(model)
%
% DESCRIPTION:
%   Launch an interactive plotting interface for the fluid model.
%
% REQUIRED PARAMETERS:
%   model - Some `ReservoirModel`-derived class with a valid fluid model.
%           This function has primarily been tested for black-oil and
%           black-oil similar fluids.
%
% OPTIONAL PARAMETERS:
%   'pressureRange'   -  An array of the pressures values to be used for
%                        pressure-dependent properties. Defaults to
%                        `max(model.minimumPressure, 0.1*barsa)` to 
%                        `min(model.maximumPressure, 600*barsa)`
% RETURNS:
%   h    - Figure handle to the plotting panel.
%

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
    opt = struct('pressureRange', [], 'regNo', 1, 'field', []);
    opt = merge_options(opt, varargin{:});
    fn = fieldnames(model.fluid);
    for propno = 1:numel(fn)
        f = fn{propno};
        lfn = model.fluid.(f);
        if iscell(lfn)
            index = min(numel(lfn), opt.regNo);
            fprintf('Multiple regions found for %s. Only plotting region %d.\n', f, index);
            model.fluid.(f) = lfn{index};
        end
    end
    if isempty(opt.pressureRange)
        p0 = max(model.minimumPressure, 0.1*barsa);
        p1 = min(model.maximumPressure, 600*barsa);
        opt.pressureRange = subdiv(p0, p1);
    else
        opt.pressureRange = reshape(opt.pressureRange, [], 1);
    end
    % Make a figure that's wider than the default.
    df = get(0, 'DefaultFigurePosition');
    fh = figure('Position', df.*[1 1 1.75 1]);
    if nargout > 0
        varargout{1} = fh;
    end
    % We want to ensure that the lines are nice and pretty even if the user
    % has messed with the defaults.
    set(fh, 'Renderer', 'painters');
    
    % Options can alter the amount of space the plot itself takes up.
    lm = 0.1;
    pw = .7;
    % Somewhat magic numbers because the gui in matlab has some magic
    % constants itself.
    plotaxis  = subplot('Position', [.75*lm, lm, pw, 1-2*lm]);
    ctrlpanel = uipanel('Position', [lm+pw, .25*lm, 1-1.25*lm-pw, 1-1.25*lm], ...
                        'Title',  'Property');
    true3ph = model.water && model.oil && model.gas;
    hasPolymer = isprop(model, 'polymer') && model.polymer;
    names = {'Viscosity',...
             'Relative permeability', ...
             'Rock compressibility', ...
             'Densities', ...
             'Capillary pressure', ...
             'Rs and Rv', ...
             '3ph-relperm: Water', ...
             '3ph-relperm: Oil', ...
             '3ph-relperm: Gas', ...
             'Polymer adsorption', ...
             'Polymer Water viscosity multiplier'...
             };
    functions = {@(name) plotViscosity(model, name), ...
                 @(name) plotRelperm(model, name),...
                 @(name) plotPVMult(model, name),...
                 @(name) plotDensity(model, name), ...
                 @(name) plotCapillaryPressure(model, name), ...
                 @(name) plotMaxR(model, name), ...
                 @(name) plot3phRelPerm(model, name, 1), ...
                 @(name) plot3phRelPerm(model, name, 2), ...
                 @(name) plot3phRelPerm(model, name, 3), ...
                 @(name) plotPolyAdsorption(model, name), ...
                 @(name) plotPolyViscMult(model, name)...
                };
    
    disgas = isprop(model, 'disgas') && model.disgas;
    vapoil = isprop(model, 'vapoil') && model.vapoil;
    active = [true; ...
              true; ...
              isfield(model.fluid, 'pvMultR'); ...
              true; ...
              isfield(model.fluid, 'pcOW') || isfield(model.fluid, 'pcOG'); ...
              disgas || vapoil; ...
              true3ph;
              true3ph;
              true3ph;
              hasPolymer && isfield(model.fluid, 'ads');
              hasPolymer && isfield(model.fluid, 'muWMult')];
              
              
    names = names(active);
    functions = functions(active);
    
    phases = {'Water', 'Oil', 'Gas'};
    for it = 1:numel(phases)
        ph = phases{it};
        if model.(lower(ph))
            names = [names, [ph, ' viscosity'], ...
                            [ph, ' b-factor']];                            %#ok<AGROW>
            functions = [functions, ...
               {@(name) plot2D_property(model, {['mu', ph(1)]}, name), ...
                @(name) plot2D_property(model, {['b', ph(1)]}, name)}];          %#ok<AGROW>
        end
    end
    
    propsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'listbox',...
              'String', names, 'Callback', @drawPlot, ...
              'Position',[0 0 1 1]);
          
    if ~isempty(opt.field)
        v = find(strcmpi(opt.field, names));
        if isempty(v), v = 1; end
        set(propsel, 'Value', v);
    end

    function drawPlot(src, event)                                          %#ok<*INUSD>
        axis(plotaxis);
        colorbar off;
        axis normal
        view(0, 90);
        ix = get(propsel, 'Value');
        fn = functions{ix};
        name = names{ix};
        
        fn(name);
    end

    drawPlot([], []);


    function plot2D_property(model, fields, plottitle)
        axis on
        f = model.fluid;
        p = opt.pressureRange;
        nval = numel(p);
        s = subdiv(0, 1, nval);
        arg = {};

        rsMax = 0;
        rvMax = 0;
        if disgas
            rsMax = f.rsSat(p, arg{:});
        end
        if vapoil
            rvMax = f.rvSat(p, arg{:});
        end

        legflag = false(size(fields));
        legh = zeros(size(fields));
        ctr = 0;
        yl = '';
        cla;
        hold on
        nf = numel(fields);
        colors = lines(nf);
        
        for i = 1:nf
            fn = fields{i};
            legflag(i) = true;
            ctr = ctr + 1;
            switch(lower(fn))
                case {'krw', 'krg', 'krog', 'kro', 'krow'}
                    data = f.(fn)(s, arg{:});
                    x = s;
                    xl = 'Saturation';
                case {'pcow', 'pcog'}
                    data = f.(fn)(s, arg{:});
                    x = s;
                    data = data/barsa;
                    xl = 'Saturation';
                    yl = 'Capillary pressure [bar]';
                case {'bo'}
                    [x, data, ok] = evalSat(model, f, fn, p, rsMax, rvMax, arg{:});
                    x = x/barsa;
                    yl = 'Shrinkage factor b(p) = 1/B(p)';
                    xl = 'Pressure [bar]';
                case {'bg', 'bw'}
                    [x, data, ok] = evalSat(model, f, fn, p, rsMax, rvMax, arg{:});
                    x = x/barsa;
                    yl = 'Expansion factor b(p) = 1/B(p)';
                    xl = 'Pressure [bar]';
                case {'rhoo', 'rhog', 'rhow'}
                    phLetter = fn(end);
                    bsub = ['b', phLetter];
                    rho = f.(['rho', phLetter, 'S']);
                    [x, b, ok, r] = evalSat(model, f, bsub, p, rsMax, rvMax, arg{:});
                    if strcmpi(phLetter, 'o') && disgas
                        % Account for solution gas
                        rhoO = rho;
                        rhoG = f.rhoGS;
                        data = b.*(r.*rhoG + rhoO);
                    elseif strcmpi(phLetter, 'g') && vapoil
                        % Account for vaporized oil
                        rhoG = rho;
                        rhoO = f.rhoOS;
                        data = b.*(r.*rhoO + rhoG);
                    else
                        data = b*rho;
                    end
                    x = x/barsa;
                    xl = 'Pressure [bar]';
                    yl = 'Density [kg/m^3]';
                case {'muw', 'muo', 'mug'}
                    [x, data, ok] = evalSat(model, f, fn, p, rsMax, rvMax, arg{:});
                    data = data/(centi*poise);
                    x = x/barsa;
                    xl = 'Pressure [bar]';
                    yl = 'Viscosity [cP]';
                case {'rssat', 'rvsat'}
                    data = f.(fn)(p, arg{:});
                    x = p/barsa;
                    xl = 'Pressure [bar]';
                case {'pvmultr'}
                    data = f.(fn)(p);
                    x = p/barsa;
                    xl = 'Pressure [bar]';
                    yl = 'Pore volume multiplier';
                case 'ads'
                    x = (0:0.01:f.cmax)';
                    data = f.(fn)(x, arg{:});
                    yl = 'Adsorption';
                    xl = 'Polymer concentration';
                case 'muwmult'
                    x = (0:0.01:f.cmax)';
                    data = f.(fn)(x, arg{:});
                    yl = 'Viscosity multiplier';
                    xl = 'Polymer concentration';
            end
            if size(data, 2) > 1
                for j = 1:size(data, 2)
                    o = ok(:, j);
                    % Ok are saturated lines, draw as thick lines
                    h = plot(x(o, j), data(o, j), 'linewidth', 2, 'color', colors(i, :));
                    % Not ok are undersaturated, draw as thin lines
                    plot(x(~o, j), data(~o, j), '--', 'linewidth', 1, 'color', colors(i, :));
                end
            else
                h = plot(x, data, 'linewidth', 2, 'color', colors(i, :));
            end
            legh(i) = h(1);
        end
        
        grid on
        legend(legh, fields(legflag),'Location','Best')
        xlabel(xl)
        ylabel(yl);
        if nargin == 3
            title(plottitle)
        end
        axis tight
    end

    function plotViscosity(model, name)
        vnames = {'muW', 'muO', 'muG'};
        act = model.getActivePhases();
        plot2D_property(model, vnames(act));
        title('Viscosity');
    end

    function plotRelperm(model, name)
        krnames = {};
        if model.water
            krnames = [krnames, 'krW'];
        end
        if model.oil
            if model.water && model.gas
                krnames = [krnames, 'krOW', 'krOG'];
            else
                if isfield(model.fluid, 'krO')
                    krnames = [krnames, 'krO'];
                else
                    krnames = [krnames, 'krOW'];
                end
            end
        end
        if model.gas
            krnames = [krnames, 'krG'];
        end
        plot2D_property(model, krnames);
        title(name);
    end

    function plotDensity(model, name)
        rnames = {'rhoW', 'rhoO', 'rhoG'};
        act = model.getActivePhases();
        plot2D_property(model, rnames(act));
        title(name);
    end

    function plotCapillaryPressure(model, name)
        cnames = {};
        fld = {'pcOW', 'pcOG'};
        for i = 1:numel(fld)
            if isfield(model.fluid, fld{i})
                cnames = [cnames, fld{i}];                                 %#ok<AGROW>
            end
        end

        plot2D_property(model, cnames);
        title('Capillary pressure');
    end

    function plotMaxR(model, name)
        rnames = {};
        if disgas
            rnames = [rnames, 'rsSat'];
        end
        if vapoil
            rnames = [rnames, 'rvSat'];
        end
        plot2D_property(model, rnames);
        title(name);
    end

    function plotPVMult(model, name)
        plot2D_property(model, {'pvMultR'});
        title(name);
    end
    
    function plot3phRelPerm(model, name, ix)                               %#ok<INUSL>
        cla;
        ns = 50;
        s = subdiv(0, 1, ns);
        [x, y] = meshgrid(s);
        [krW, krO, krG] = ...
            model.relPermWOG(x(:), 1-x(:)-y(:), y(:), model.fluid);
        
        % If sW and sG sum up to more than unity, pad with NaN.
        unphys = x + y > 1;
        krW(unphys) = nan;
        krO(unphys) = nan;
        krG(unphys) = nan;
        
        [mapx, mapy] = ternaryAxis('names', {'S_w', 'S_o', 'S_g'});
        xx = mapx(x, 1-x-y, y);
        yy = mapy(x, 1-x-y, y);
        if ix == 1
            contourf(xx, yy, reshape(krW, size(xx)),11)
            title('Water relative permeability')
        elseif ix == 2
            contourf(xx, yy, reshape(krO, size(xx)),11)
            title('Oil relative permeability')
        else
            contourf(xx, yy, reshape(krG, size(xx)),11)
            title('Gas relative permeability')
        end
        axis equal

        xlim([0, 1]);
        ylim([0, 1]);
        caxis([0, 1]);
        xlabel('S_w')
        ylabel('S_g')
        colorbar('horizontal')
        legend off
    end

    function x = subdiv(start, stop, n)
        if nargin < 3
            n = 100;
        end
        dx = (stop - start)/(n-1);
        x = (start:dx:stop)';
    end

    function [x, y, ok, rs_g] = evalSat(model, f, fn, x, rsMax, rvMax, varargin)
        ok = true(size(x));
        rs_g = nan;
        if checkBO(model)
            if any(strcmpi(fn, {'muo', 'bo'})) && disgas
                mrs = max(rsMax);
                rs = 0:mrs/20:mrs;
                [x, rs_g] = meshgrid(x, rs);
                rssat = zeros(size(x));
                for i = 1:size(x, 1)
                    rssat(i, :) = f.rsSat(x(i, :), varargin{:});
                end

                saturated = rs_g >= rssat;
                rs_g(saturated) = rssat(saturated);
                ok = saturated';
                y = f.(fn)(x, rs_g, saturated, varargin{:})';
                x = x';
                rs_g = rs_g';
            elseif any(strcmpi(fn, {'mug', 'bg'})) && vapoil
                mrs = max(rvMax);
                rv = 0:mrs/10:mrs;
                [x, rs_g] = meshgrid(x, rv);
                rvsat = zeros(size(x));
                for i = 1:size(x, 1)
                    rvsat(i, :) = f.rvSat(x(i, :));
                end

                saturated = rs_g >= rvsat;
                rs_g(saturated) = rvsat(saturated);
                ok = saturated';
                y = f.(fn)(x, rs_g, saturated, varargin{:})';
                x = x';
                rs_g = rs_g';
            else
                y = f.(fn)(x, varargin{:});
                rs_g = nan;
            end
        else
            y = f.(fn)(x, varargin{:});
        end
    end

    function plotPolyAdsorption(model, name)
        plotStuff(model, {'ads'});
        title(name);
    end

    function plotPolyViscMult(model, name)
        plotStuff(model, {'muWMult'});
        title(name);
    end

    function ind = checkBO(model)
        ind = isa(model, 'ThreePhaseBlackOilModel') &&...
               (disgas || vapoil);
    end
end