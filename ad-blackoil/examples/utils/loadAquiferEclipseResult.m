function [statesEcl, G] = loadEclipseResult(dir, fnroot)

    fn = fullfile(dir, [fnroot, '.INIT']);
    init = readEclipseOutputFileUnFmt(fn);
    fn = fullfile(dir, [fnroot, '.EGRID']);
    grid = readEclipseOutputFileUnFmt(fn);
    [eclG, eclrock, N, T] = initGridFromEclipseOutput(init, grid, 'outputSimGrid', ...
                                                      true);
    G = eclG{1};

    statesEcl = {};
    for i = 1 : 21
        fnname = sprintf('%s.X%04.0f', fnroot, i);
        fn = fullfile(dir, fnname);
        stateEclFmt = readEclipseOutputFileUnFmt(fn);
        stateEcl.pressure = (stateEclFmt.PRESSURE.values);
        sw = stateEclFmt.SWAT.values;
        stateEcl.s = [sw, (1 - sw)];
        statesEcl{i} = stateEcl;
    end
end
