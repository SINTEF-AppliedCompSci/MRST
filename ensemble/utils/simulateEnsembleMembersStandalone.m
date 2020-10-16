function simulateEnsembleMembersStandalone(fileName, range)
    data = load(fileName);
    data.ensemble.simulateEnsembleMembersCore(range);
end