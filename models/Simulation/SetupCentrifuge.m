function model = SetupCentrifuge(model)
    G          = model.grid.G;
    rock       = model.rock;
    fluid      = model.fluid;
    g          = model.gravity;
    simulation = model.simulation;
    nCells     = G.cells.num;
    
    if (isfield(simulation,'bCells'))
        coreLength = model.experiment.geometry.length.value;
    else
        coreLength = model.experiment.geometry.length.value * (1 - 2 / nCells);
    end
    
    centRad = model.experiment.schedule.centRad.value; 
    highgrav = g * ((centRad + coreLength/2)/centRad);
    lowgrav  = g * ((centRad - coreLength/2)/centRad);
    dgz  = linspace(lowgrav,highgrav,nCells - 1)';
    dgz(1) = 0; dgz(end) = 0; 
    
    gravity('on','x', g);
    model.twoPhaseOilWaterModel = TwoPhaseOilWaterModel(G, rock, fluid);
    model.twoPhaseOilWaterModel.operators.gdz = zeros(nCells-1,1);
    Grad = model.twoPhaseOilWaterModel.operators.Grad(model.twoPhaseOilWaterModel.G.cells.centroids);
    model.twoPhaseOilWaterModel.operators.gdz(1:end) = Grad(:,1) .* dgz;
end