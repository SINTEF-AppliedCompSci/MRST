function levels = makeSpatialUpscalingLevels(setup, varargin)
    opt = struct('numLevels', 3);
    [opt, extra] = merge_options(opt, varargin{:});
    model = setup.model;
    partitions0 = makeNestedPartitions2(model, 2*opt.numLevels-1, extra{:});
    partitions0 = partitions0(1:2:end);
    partitions  = partitions0;
    for i = 2:opt.numLevels
        level = repmat(i, setup.model.G.cells.num, 1);
        level(vertcat(setup.schedule.control(1).W.cells)) = 1;
        partitions{i} = makeHybridPartition(partitions0, level, inf, []);
    end
    partitions = flipud(partitions);
    
    levels = cell(opt.numLevels, 1);
    for i = 1:opt.numLevels
        levels{i} = SpatialUpscalingLevel(i, partitions{i});
    end
end