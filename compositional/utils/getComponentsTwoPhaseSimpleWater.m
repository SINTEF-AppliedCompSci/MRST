function components = getComponentsTwoPhaseSimpleWater(model, rho, sW, xM, yM)
    if model.water
        components = cellfun(@(x, y) {[], rho{2}.*x, rho{3}.*y}, xM, yM, 'UniformOutput', false);
        components = [components, {{rho{1}.*sW, [], []}}];
    else
        components = cellfun(@(x, y) {rho{1}.*x, rho{2}.*y}, xM, yM, 'UniformOutput', false);
    end
end
