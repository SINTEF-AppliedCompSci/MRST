function [totalGasVolume, H2InWater, vH2InGas] = computeGasVolumes(model, states, steps)
    % Get phase names and find indices for liquid and gas
    names = model.getPhaseNames();
    nph = model.getNumberOfPhases;               
    for ph = 1:nph    
        switch names(ph)
            case 'O'
                L_ix = ph;
            case 'G'
                V_ix = ph;
        end
    end
    % Preallocate arrays for efficiency
    H2InWater = cell(length(steps), 1);
    vH2InGas = cell(length(steps), 1);
    totalGasVolume = cell(length(steps), 1);

    % Loop through each step to compute gas volumes
    for i = 1:length(steps) 
        step = steps(i);
        H2InWater{i} = states{step}.PVTProps.ShrinkageFactors{L_ix} .* ...
                          states{step}.s(:, V_ix) .* ...
                          model.G.cells.volumes .* ...
                          states{step}.rs;
        vH2InGas{i} = states{step}.PVTProps.ShrinkageFactors{V_ix} .* ...
                         model.G.cells.volumes .*states{step}.s(:, V_ix);
        totalGasVolume{i} = H2InWater{i} + vH2InGas{i}; % Total gas volume in the reservoir
    end
end
