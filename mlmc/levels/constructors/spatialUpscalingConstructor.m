function setup = spatialUpscalingConstructor(setup, partition, varargin)

    opt = struct('cartBox', false);
    opt = merge_options(opt, varargin{:});

    nc = setup.model.G.cells.num;

    model = upscaleModelTPFA(setup.model, partition, 'validatePartition', false);
    if opt.cartBox
        model.G.cartDims = ones(1, model.G.physDims)*sqrt(model.G.cells.num, 1);
    end
    % Map well cells
    schedule = setup.schedule;
    for i = 1:numel(schedule.control)
        W = setup.schedule.control(i).W;
        if ~isempty(W)
            for j = 1:numel(W)
               W(j).cells = model.G.partition(W(j).cells);
            end
            schedule.control(i).W = W;
        end
    end
    % Map state
    fnames = fieldnames(setup.state0);
    donor  = setup.state0;
    state0 = struct();
    for i = 1:numel(fnames)
        n = fnames{i};
        if isa(donor.(n), 'double') && size(donor.(n), 1) == nc
            state0.(n) = zeros(model.G.cells.num, size(donor.(n),2));
            state0.(n)(model.G.partition,:) = donor.(n);
        end
    end

    setup.state0 = state0;
    setup.model = model;
    setup.schedule = schedule;

end