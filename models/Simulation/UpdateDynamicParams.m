function dynamic = UpdateDynamicParams(model)
    bc          = model.bc;
    static      = model.static;
    dynamic     = model.dynamic;
    params      = dynamic.params;  
    Swi         = model.experiment.rock.Swi.value;
    simulation  = model.simulation;
    if (isfield(simulation,'bCells'))
        pv      = sum(model.rock.pv(2:end-1));
    else
        pv      = sum(model.rock.pv(1:end));
    end
    process     = model.experiment.process;
    processName = process.name;
    processType = lower(process.type);    
    bcl         = bc.left;            
    G           = model.grid.G;
    coreLength = model.experiment.geometry.length.value;
    x = [0;G.cells.centroids(2:end-1,1) - 2 * G.cells.centroids(1,1);coreLength];


    for i = 1:length(params.scheduleSteps)
        state  = model.state{i, 1};
        waterSatFull = [interp1(G.cells.centroids(2:3,1),state.s(2:3,1),0,'linear','extrap');...
            state.s(2:end-1,1);interp1(G.cells.centroids(end-2:end-1,1) - ...
            G.cells.centroids(1,1),state.s(end-2:end-1,1),coreLength,'linear','extrap')];
        if (isfield(simulation,'bCells'))
            Sw = state.s(2:end-1,1);
        else
            Sw = state.s(:,1);
        end               
        
        % injected fluids
        water_injection_rate = params.qinj(:,1);
        oil_injection_rate = params.qinj(:,2); 
        water_injection_volume = params.Qinj(:,1);
        oil_injection_volume = params.Qinj(:,2);            
        
        if(processType == "cent")
            water_injection_rate = [water_injection_rate; 0]; injection_rate(:,1) = water_injection_rate;
            oil_injection_rate = [oil_injection_rate; 0]; injection_rate(:,2) = oil_injection_rate;
            water_injection_volume = [water_injection_volume; 0]; injection_volume(:,1) = water_injection_volume;
            oil_injection_volume = [oil_injection_volume; 0]; injection_volume(:,2) = oil_injection_volume;
        else                       
            water_injection_rate = [water_injection_rate; bcl.qw_inj]; injection_rate(:,1) = water_injection_rate;
            oil_injection_rate = [oil_injection_rate; bcl.qo_inj]; injection_rate(:,2) = oil_injection_rate;        
            water_injection_volume = [water_injection_volume; water_injection_volume(end) + water_injection_rate(end) * params.scheduleSteps(i)];            
            oil_injection_volume = [oil_injection_volume; oil_injection_volume(end) + oil_injection_rate(end) * params.scheduleSteps(i)];
            injection_volume(:,1) = water_injection_volume; injection_volume(:,2) = oil_injection_volume;
        end
        params.qinj = injection_rate; 
        params.Qinj = injection_volume;
        clear injection_rate injection_volume

        PVI = params.PVI;
        params.PVI = [PVI; (water_injection_volume(end) + oil_injection_volume(end)) / pv];

        % produced fluids     
        water_production_rate = params.qp_net(:,1);
        oil_production_rate = params.qp_net(:,2); 
        water_production_volume = params.Qp_net(:,1);
        oil_production_volume = params.Qp_net(:,2);
        F = griddedInterpolant(x,waterSatFull);
        fun = @(t) F(t);
        q = integral(fun, x(1), x(end))/max(x);
        if(strcmpi(processName,"imbibition"))
            water_production_volume = [water_production_volume; 0]; net_production_volume(:,1) = water_production_volume;
            water_production_rate = [water_production_rate; water_production_rate(end) / params.scheduleSteps(i)]; net_production_rate(:,1) = water_production_rate;
            oil_production_volume = [oil_production_volume; pv*((1-Swi)-(1-q))]; net_production_volume(:,2) = oil_production_volume;
            oil_production_rate = [oil_production_rate; oil_production_rate(end) / params.scheduleSteps(i)]; net_production_rate(:,2) = oil_production_rate;
        elseif(strcmpi(processName,"drainage"))
            water_production_volume = [water_production_volume; pv*(Swi-(q))]; net_production_volume(:,1) = water_production_volume;
            water_production_rate = [water_production_rate; water_production_rate(end) / params.scheduleSteps(i)]; net_production_rate(:,1) = water_production_rate;
            oil_production_volume = [oil_production_volume; 0]; net_production_volume(:,2) = oil_production_volume;
            oil_production_rate = [oil_production_rate; oil_production_rate(end) / params.scheduleSteps(i)]; net_production_rate(:,2) = oil_production_rate;
        end 

        % Net production
        params.qp_net = net_production_rate;
        params.Qp_net = net_production_volume;
        
        % Cumulative production
        params.qprod = net_production_rate + [water_injection_rate, oil_injection_rate];
        params.Qprod = net_production_volume  + [water_injection_volume, oil_injection_volume];
        clear net_production_rate net_production_volume
        
        pDiff = params.pDiff;
        params.pDiff = [pDiff;state.pressure(static.gauge.left)-...
                              state.pressure(static.gauge.right)];            
        params.gradp = gradient(pDiff);

        SwAvg = params.SwAvg;
        params.SwAvg = [SwAvg;mean(Sw)];
        if min(Sw) < params.Sw_min
            params.Sw_min = min(Sw);
        end
        if max(Sw) > params.Sw_max
            params.Sw_max = max(Sw);
        end

        dynamic.params = params;
        model.dynamic = dynamic;        
    end

end