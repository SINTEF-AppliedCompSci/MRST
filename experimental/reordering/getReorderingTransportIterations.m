function iterations = getReorderingTransportIterations(reports)

    iterations = cell(numel(reports),1);
    [iterations{:}] = deal(0);
    for i = 1:numel(reports)
        srOuter = reports{i}.StepReports;
        for j = 1:numel(srOuter)
            sr = srOuter{j}.NonlinearReport{1}.TransportSolver.StepReports;
            for k = 1:numel(sr)
                iterations{i} = iterations{i} + sr{k}.NonlinearReport{1}.Iterations;
            end
        end
    end

end