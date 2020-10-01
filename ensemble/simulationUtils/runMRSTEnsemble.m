function runMRSTEnsemble(ensemble, vargin)
    % Run simulation for all the members in the MRSTEnsemble ensemble 
    
    opts = struct('clearPackedSimulationOutput', 'true');
    [opts, extra] = merge_options(opts, vargin{:});
    
    for seed=1:ensemble.num
        ensemble.simulateEnsembleMember(seed, opts);
    end
    
end

