function fluidPlotPanelAD(model)
    % Make a figure that's wider than the default.
    df = get(0, 'DefaultFigurePosition');
    fh = figure('Position', df.*[1 1 1.75 1]);
    % We want to ensure that the lines are nice and pretty even if the user
    % has messed with the defaults.
    set(fh, 'Renderer', 'painters');
    
    % Options can alter the amount of space the plot itself takes up.
    lm = 0.1;
    pw = .7;
    % Somewhat magic numbers because the gui in matlab has some magic
    % constants itself.
    plotaxis  = subplot('Position', [.75*lm, lm, pw, 1-2*lm]);
    ctrlpanel = uipanel('Position', [lm+pw, .25*lm, 1-1.25*lm-pw,  1-1.25*lm], ...
                        'Title',  'Property');
    
    names = {'Viscosity',...
             'Relperm', ...
             'Rock compressibility', ...
             'Compressibility'};
    if isfield(model.fluid, 'pcOW') || isfield(model.fluid, 'pcOG')
        names{end + 1} = 'Capillary pressure';
    end
    if checkBO(model)
        names{end + 1} = 'Max dissolution';
    end
    propsel = uicontrol('Units', 'normalized', 'Parent', ctrlpanel,...
              'Style', 'listbox',...
              'String', names, 'Callback', @drawPlot, ...
              'Position',[0 0 1 1]);
    
          

    function drawPlot(src, event)
        axis(plotaxis);
        name = names{get(propsel, 'Value')};
        switch lower(name)
            case 'viscosity'
                plotViscosity(model, name)
            case 'relperm'
                plotRelperm(model, name)
            case 'capillary pressure'
                plotCapillaryPressure(model, name)
            case 'compressibility'
                plotCompressibility(model, name)
            case 'max dissolution'
                plotMaxR(model, name);
            case 'rock compressibility'
                plotPVMult(model, name);
            otherwise
                disp(name)
        end
    end

    drawPlot([], []);
end

function plotStuff(model, fields)
    f = model.fluid;
    p = subdiv(0, 1000*barsa);
    s = subdiv(0, 1);

    rs = 0;
    rv = 0;
    n = sum(model.getActivePhases);

    legflag = false(size(fields));
    data = nan(numel(p), n);
    
    ctr = 0;
    
    for i = 1:numel(fields)
        fn = fields{i};
        if ~isfield(f, fn)
            continue
        end
        legflag(i) = true;
        ctr = ctr + 1;
        switch(lower(fn))
            
            case {'krw', 'krg', 'krog', 'kro', 'krow'}
                data(:, ctr) = f.(fn)(s);
                x = s;
                xl = 'Saturation';
            case {'pcow', 'pcog'}
                data(:, ctr) = f.(fn)(s);
                x = s;
                xl = 'Saturation';
            case {'bo', 'bg', 'bw'}
                data(:, ctr) = evalSat(model, f, fn, p, rs, rv);
                % Subtract from lowest value to avoid scaling issues with
                % graph
                data(:, ctr) = data(:, ctr) - data(1, ctr);
                x = p/barsa;
                xl = 'Pressure (barsa)';
            case {'muw', 'muo', 'mug'}
                data(:, ctr) = evalSat(model, f, fn, p, rs, rv);
                x = p/barsa;
                xl = 'Pressure (barsa)';
            case {'rssat', 'rvsat'}
                data(:, ctr) = f.(fn)(p);
                x = p/barsa;
                xl = 'Pressure (barsa)';
            case {'pvmultr'}
                data(:, ctr) = f.(fn)(p);
                x = p/barsa;
                xl = 'Pressure (barsa)';
        end
    end
    plot(x, data, 'linewidth', 2)
    grid on
    legend(fields(legflag))
    xlabel(xl)
end

function plotViscosity(model, name)
    plotStuff(model, {'muW', 'muO', 'muG'});
    if checkBO(model)
        title([name, ' (saturated curves)']);
    else
        title(name);
    end
end

function plotRelperm(model, name)
    plotStuff(model, {'krW', 'krOW', 'krOG', 'krO', 'krG'});
    title(name);
end

function plotCompressibility(model, name)
    plotStuff(model, {'bW', 'bO', 'bG'});
    if checkBO(model)
        title([name, ' (saturated curves)']);
    else
        title(name);
    end
end

function plotCapillaryPressure(model, name)
    plotStuff(model, {'pcOW', 'pcOG'});
    title('Capillary pressure');
end

function plotMaxR(model, name)
    plotStuff(model, {'rsSat', 'rvSat'});
    title(name);
end

function plotPVMult(model, name)
    plotStuff(model, {'pvMultR'});
    title(name);
end

function x = subdiv(start, stop)
    dx = (stop - start)/250;
    x = (start:dx:stop)';
end

function y = evalSat(model, f, fn, x, rs, rv)
    if checkBO(model)
        if any(strcmpi(fn, {'muo', 'bo'})) && model.disgas
            y = f.(fn)(x, rs, true);
        elseif any(strcmpi(fn, {'mug', 'bg'})) && model.vapoil
            y = f.(fn)(x, rv, true);
        else
            y = f.(fn)(x);
        end
    else
        y = f.(fn)(x);
    end
end

function ind = checkBO(model)
    ind = isa(model, 'ThreePhaseBlackOilModel') &&...
           (model.disgas || model.vapoil);
end
