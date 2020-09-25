function runMRSTEnsemble(ensemble, vargin)
    % Run simulation for all the members in the MRSTEnsemble ensemble 
    
    opts = struct('clearPackedSimulationOutput', 'true');
    [opts, extra] = merge_options(opts, vargin{:});
    
    for i=1:ensemble.num
        ensemble.runSimulation(i, opts);
    end
    
end

