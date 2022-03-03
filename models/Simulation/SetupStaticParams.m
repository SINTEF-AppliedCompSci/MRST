function static = SetupStaticParams(model)
    if(isfield(model.simulation,'gaugeOff'))
        gaugeOff = model.simulation.gaugeOff.value;
        nCells   = model.simulation.nCells.value;
        length   = model.experiment.geometry.length.value;        
        dx = length / nCells;
        gauge.left = round(gaugeOff / dx) + 1;
        gauge.right = nCells - gauge.left + 1;       
    else
        gauge.left  = 1;
        gauge.right = model.grid.G.cells.num;
    end    
    static.gauge = gauge;    
    static.dt  = model.simulation.timeStep.value;   
    %---------------------------------------
    fig.title    = 'Visualization';
    if(isfield(model.plot,'style'))
        if(strcmp(model.plot.style.inputStyle,'normal')), fig.style = 'normal'; end
        if(strcmp(model.plot.style.inputStyle,'docked')), fig.style = 'docked'; end
    end    
    fig.tag      = fig.title;
    fig.subtitle = {'Flooding Experiment',...
                    'Injection Rates',...
                    'Saturation Front',...
                    'Pressure Differential',...
                    'Average Water Saturation',...
                    'Production'};                
    static.fig = fig;
    %---------------------------------------
end