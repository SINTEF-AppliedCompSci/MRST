function dynamic = UpdateDynamicParams_usingrates(model)
    static      = model.static;
    dynamic     = model.dynamic;
    params      = dynamic.params;  
    simulation  = model.simulation;
    process     = model.experiment.process;
    processName = process.name;
    processType = lower(process.type);    
    cumScheduleSteps = params.cumScheduleSteps;
    if( strcmpi(processType,"cent") && (strcmpi(processName,"imbibition")) ) 
        if (isfield(simulation,'bCells'))
            faceIn = model.grid.G.cells.num; faceOut = 2;
            pv          = sum(model.rock.pv(2:end-1));
        else
            faceIn = model.grid.G.cells.num + 1; faceOut = 1;
            pv          = sum(model.rock.pv);
        end
    else
        if (isfield(simulation,'bCells'))
            faceIn = 2; faceOut = model.grid.G.cells.num;
            pv          = sum(model.rock.pv(2:end-1));
        else
            faceIn = 1; faceOut = model.grid.G.cells.num + 1;
            pv          = sum(model.rock.pv);
        end
    end


    for i = 1:length(params.scheduleSteps)
        state  = model.state{i, 1};
        if (isfield(simulation,'bCells'))
            Sw = state.s(2:end-1,1);
        else
            Sw = state.s(:,1);
        end 

        % injected fluids
        water_injection_rate = params.qinj(:,1);
        oil_injection_rate = params.qinj(:,2); 
        
        water_injection_rate = [water_injection_rate; state.flux(faceIn,1)]; injection_rate(:,1) = water_injection_rate;
        oil_injection_rate = [oil_injection_rate; state.flux(faceIn,2)]; injection_rate(:,2) = oil_injection_rate;        
        params.qinj = injection_rate; 

        % produced fluids     
        water_production_rate = params.qprod(:,1);
        oil_production_rate = params.qprod(:,2); 
        water_production_rate = [water_production_rate; state.flux(faceOut, 1)]; production_rate(:,1) = water_production_rate;
        oil_production_rate = [oil_production_rate; state.flux(faceOut, 2)]; production_rate(:,2) = oil_production_rate;

        % Net production
        params.qp_net = production_rate - injection_rate;
        
        
%       % Cumulative production
        params.qprod = production_rate ;
        
        pDiff = params.pDiff;
        params.pDiff = [pDiff;state.pressure(static.gauge.left)-...
                              state.pressure(static.gauge.right)];            
        params.gradp = gradient(pDiff);

        SwAvg = params.SwAvg;
        params.SwAvg = [SwAvg;mean(Sw)];
        
        clear injection_rate production_rate
    
    end
        
    water_injection_volume = cumtrapz(cumScheduleSteps, params.qinj(:,1));            
    oil_injection_volume = cumtrapz(cumScheduleSteps, params.qinj(:,2));
    injection_volume(:,1) = water_injection_volume; injection_volume(:,2) = oil_injection_volume;
    params.Qinj = injection_volume;
    
    params.PVI = (water_injection_volume + oil_injection_volume) ./ pv;
    
    water_production_volume = cumtrapz(cumScheduleSteps, params.qprod(:,1)); production_volume(:,1) = water_production_volume;           
    oil_production_volume = cumtrapz(cumScheduleSteps, params.qprod(:,2)); production_volume(:,2) = oil_production_volume;
    params.Qp_net = production_volume - injection_volume ;
    params.Qprod = production_volume;
        
    dynamic.params = params;
    model.dynamic = dynamic;        
end