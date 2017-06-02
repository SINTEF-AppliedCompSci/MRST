function example4(varargin)

    opt.recompute     = false;
    opt.savedir       = 'results/example4';
    opt.produce_plots = true;
    opt.subtrap_file  = 'utsira_subtrap_function.mat';
    opt.inj_steps     = 10;
    opt.inj_years     = 50;
    opt.migr_steps    = 60;
    opt.migr_years    = 3000;
    opt.dis_rate      = 8.6924e-11;
    opt.dis_max       = 0.02;
    
    opt = merge_options(opt, varargin{:});
    mrstModule add ad-core;
    gravity on;
    
    %% Loading grid and subscale trapping function
    [Gt, rock2D] = getUtsiraTopGrid(1, true);
    dh = [];
    if ~isempty(opt.subtrap_file)
        recenter = true;
        dh = computeToptrapsUtsira(opt.subtrap_file, Gt, recenter);
    end
    
    %% Defining injection point and schedules
    ref_rho         = getValuesSPE134891();
    totvols         = [2 10] * 1e3 * 1e6 * kilogram * opt.inj_years / ref_rho(2);
    sleipner_coords = [438516, 6471210];
    wcell           = closest_gridcell(Gt, sleipner_coords);
    Winj            = addWellVE([], Gt, rock2D, wcell , ...
                                'type'   , 'rate'     , ...
                                'val'    , 1          , ... % correct value set later
                                'radius' , 0.3        , ...
                                'comp_i'  , [0 1]);
    Wmig            = Winj; Wmig.val = 0;
    istep           = ones(opt.inj_steps, 1) * ...
                          (opt.inj_years * year / opt.inj_steps);
    mstep           = [istep(1); ...
                       ones(opt.migr_steps-1, 1) * ...
                       (opt.migr_years * year / (opt.migr_steps - 1))];
    mstep(2)        = mstep(2) - mstep(1);
    bfaces          = identifyBoundaryFaces(Gt);
    bcond           = addBC([], bfaces, ...
                            'pressure', Gt.faces.z(bfaces) * ref_rho(1) * ...
                                        norm(gravity) + 1 * atm, ...
                            'sat', [1 0]);
    for i = 1:2

       Winj.val = totvols(i) / (opt.inj_years * year);
       Winj.name = ['I', num2str(i)];

       schedules(i) = ...
           struct('control', [struct('W', Winj, 'bc', bcond), ...
                              struct('W', Wmig, 'bc', bcond)], ...
                  'step'   , struct('control', [1 * ones(size(istep)); ...
                                                2 * ones(size(mstep))], ...
                                    'val', [istep; mstep]));%#ok
    end
    
    injcase = {'TwoMT', 'TenMT'};
    
    %% Computing the different cases, if not already done,
    %  or if recompute is requested
    
    for dis = {'nodiss', 'instdiss', 'ratediss'} 
        for i = 1:2 
            result_dir = fullfile(opt.savedir, dis{1}, injcase{i}, 'report'); 
            if ~computed(result_dir) || opt.recompute
               
                % Compute case here (setting 'muCO2' to zero below activates
                % variable CO2 viscosity computation).
                runSimulation(Gt, rock2D, schedules(i)      ,                 ...
                              'dh'          , dh            ,                 ...
                              'report_dir'  , result_dir    ,                 ...
                              'do_plot'     , false         ,                 ...
                              'dissolution' , disval(dis    , opt.dis_rate) , ...
                              'dis_max'     , opt.dis_max   ,                 ...
                              'refRhoCO2'   , ref_rho(2)    ,                 ...
                              'refRhoBrine' , ref_rho(1)    ,                 ...
                              'muCO2'       , 0             ,                 ... 
                              'muBrine'     , 8e-4 * Pascal * second);
            end
        end
    end
    
    %% Plot result, if requested
    if opt.produce_plots
        for dis = {'nodiss', 'instdiss', 'ratediss'}
            for i = 1:2
                close all; % prevent too many windows accumulating
                
                fprintf('Current case: %s - %s\n', dis{1}, injcase{i});
                
                basedir = fullfile(opt.savedir, dis{1});
                
                [Gt, reports] = load_data(fullfile(basedir, injcase{i}, ...
                                                   'report'), 1:(opt.inj_steps ...
                                                                 + opt.migr_steps));

                tsteps = [opt.inj_steps, 31, 70];
                
                selectedResultsMultiplot(Gt, reports, tsteps, ...
                                         'background' , 'totalCO2' ,  ...
                                         'plot_traps' , true);
    
                fprintf('Press ''enter'' to advance to the next simulation case.\n\n');
                pause; 
            end
        end
    end
end

% ----------------------------------------------------------------------------

function [Gt, stepdata] = load_data(dir, tsteps)
    
    % loading grid
    tmp = load(fullfile(dir, 'simulation_grid'), 'Gt');
    Gt = tmp.Gt;
    
    % loading timesteps
    for i = 1:numel(tsteps)
       tmp = load(fullfile(dir, sprintf('report_%i', tsteps(i))), 'report');
       stepdata(i) = tmp.report;%#ok
    end
end



% ----------------------------------------------------------------------------
function res = computed(dirname)
    
    if exist(dirname) ~= 7%#ok
        res = false; return;
    end
    
    d = dir(dirname);
    filenum = sum(~[d.isdir]);
    
    res = (filenum > 0);
end
% ----------------------------------------------------------------------------
% function res = last_step(dirname)
%     % Return the number of the last saved step encountered
%     d   = dir(fullfile(dirname, 'report_*.mat'));
%     res = max(cellfun(@(x) str2num(x(8:end-4)), {d.name}));
% end

% ----------------------------------------------------------------------------
function val = disval(dis, rate)
    if strcmpi(dis, 'nodiss')
        val = 0;
    elseif strcmpi(dis, 'instdiss')
        val = Inf;
    else
        val = rate;
    end
end
% ----------------------------------------------------------------------------
function c_ix = closest_gridcell(Gt, coords)

    dvec = bsxfun(@minus, Gt.cells.centroids, coords);
    dist = sum((dvec.^2), 2);
    [~, c_ix] = min(dist);
    
end

% ----------------------------------------------------------------------------
function bfaces = identifyBoundaryFaces(Gt)
   
    bfaces = find(any(Gt.faces.neighbors==0,2)); 

end
            