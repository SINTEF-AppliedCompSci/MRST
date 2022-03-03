function dynamic = UpdateDynamicParams_3D(model)
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

    for i = 1:length(params.scheduleSteps)
        state  = model.state{i, 1};
        if (isfield(simulation,'bCells'))
            Sw = state.s(and(not(G.inlet_mask),not(G.outlet_mask)),1);
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
        x = linspace(0,numel(Sw),numel(Sw));
        integral = cumtrapz(x,Sw);
        q = integral(end) / x(end);
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
        
        inlet_mask = G.inlet_mask;
        outlet_mask = G.outlet_mask;
        
        params.pDiff = [pDiff;mean(state.pressure(inlet_mask))-...
                              mean(state.pressure(outlet_mask))];            
        params.gradp = gradient(pDiff);

        SwAvg = params.SwAvg;
        params.SwAvg = [SwAvg;mean(Sw)];

        dynamic.params = params;
        model.dynamic = dynamic;        
    end

end