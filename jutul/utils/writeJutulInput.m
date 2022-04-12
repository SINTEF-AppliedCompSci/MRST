function pth = writeJutulInput(state0, model, schedule, name, varargin)
    opt = struct('path', tempdir(), 'printcmd', true, 'extra', {{}});
    [opt, other] = merge_options(opt, varargin{:});
    assert(isfolder(opt.path))
    pth = writeToJutul(opt.path, name, state0, model, schedule, opt.extra, other);
    if opt.printcmd
        fprintf('Case written. To simulate, call:\n\n\tusing JutulDarcy\n\tsimulate_mrst_case("%s", write_mrst = true);\n\n', pth)
    end
end

function pth = writeToJutul(folder_path, name, state0, model, schedule, extra, other)
    G = model.G;
    rock = model.rock;
    rock.poro = model.operators.pv./G.cells.volumes;
    deck = model.inputdata;
    % Make sure everything is prepped as if we were running a simulation.
    model = model.validateModel();
    grav = norm(model.gravity) > 0;
    if ~grav
        G.cells.centroids(:, 3) = 0;
    end
    schedule = model.validateSchedule(schedule);
    % Create tables etc if not present as input deck.
    % Note: There is the possibility of inconsistencies here if the model
    % was modified after construction from input file.
    if isempty(deck)
        deck = struct('PROPS', [], 'RUNSPEC', []);
    end
    deck.PROPS = generatePROPS(model, 'props', deck.PROPS, 'writeExtra', true, other{:});
    deck.RUNSPEC = generateRUNSPEC(model, deck.RUNSPEC);
    assert(~isempty(deck));
    if isfield(schedule.control, 'W')
        for i = 1:numel(schedule.control)
            ctrl = schedule.control(i);
            if isfield(ctrl.W, 'lims')
                for j = 1:numel(ctrl.W)
                    ctrl.W(j) = prepWell(ctrl.W(j), grav);
                end
                schedule.control(i) = ctrl;
            end
        end
    end
    for i = 1:numel(schedule.control)
        schedule.control(i).W = applyFunction(@(x) x, schedule.control(i).W);
    end
    schedule.control = applyFunction(@(x) x, schedule.control);
    if isprop(model, 'disgas') && model.disgas
        mass = value(model.getProp(state0, 'ComponentTotalMass')');
        state0.zg = mass(:, 3)./sum(mass(:, 2:end), 2);
    end
    fprintf('Starting write of case %s...', name)
    pth = fullfile(folder_path, [name, '.mat']);
    jutul = struct();
    jutul.T = model.operators.T;
    jutul.N = model.operators.N;
    jutul.phases = model.getActivePhases();
    jutul.name = name;
    jutul.G = G;
    jutul.rock = rock;
    jutul.state0 = state0;
    jutul.deck = deck;
    jutul.schedule = schedule;
    jutul.extra = extra;
    if isa(model, 'ThreePhaseCompositionalModel')
        [jutul.eos, jutul.mixture] = getCompositionalOutputs(model);
    end
    save(pth, '-struct', 'jutul');
    disp(' ok.')
    pth(strfind(pth,'\'))='/';
end

function w = prepWell(w, grav)
    if isfield(w.lims, 'thp')
        w.lims = rmfield(w.lims, 'thp');
    end
    if ~grav
        w.refDepth = 0;
        w.dZ = 0*w.dZ;
    end
end

function [eos_info, mixture_s] = getCompositionalOutputs(model)
    eos = model.EOSModel;
    mixture = eos.CompositionalMixture;

    bic = mixture.getBinaryInteraction();
    if ~any(bic(:))
        bic = [];
    end
    eos_info = struct('name', eos.shortname(), 'volume_shift', eos.PropertyModel.volumeShift, 'bic', bic);
    n = numel(mixture.names);
    comps = applyFunction(@(i) struct('pc', mixture.Pcrit(i), 'Tc', mixture.Tcrit(i), 'Vc', mixture.Vcrit(i), 'acf', mixture.acentricFactors(i), 'mw', mixture.molarMass(i), 'name', mixture.names(i)), 1:n);
    mixture_s = struct('components', comps, 'names', mixture.names);
end