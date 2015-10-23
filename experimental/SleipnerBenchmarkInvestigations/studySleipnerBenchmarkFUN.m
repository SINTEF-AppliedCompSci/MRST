function varargout = studySleipnerBenchmarkFUN( varargin )
% study the Sleipner benchmark which comes from Singh et al 2010 (SPE
% paper).

% SYNOPSIS:
%
% studySleipnerBenchmarkFUN() % TODO: implement such that it lists all
%                               possible case options?
% studySleipnerBenchmarkFUN('gridname','IEAGHGmodel')
% studySleipnerBenchmarkFUN('refineLevel', 2)
% studySleipnerBenchmarkFUN('ratecase', 'original')
% studySleipnerBenchmarkFUN('modifyPerm', true)
% studySleipnerBenchmarkFUN('modifyPerm', true, 'perm_fac', 2.75)

% [ opt, var, model, schedule, initState ] = studySleipnerBenchmarkFUN()
% [ opt, var, model, schedule, initState, wellSols, states, sim_report ] = studySleipnerBenchmarkFUN()


% PARAMETERS:


% OPTIONS:


% RETURNS:
%   opt         (optional)
%   var         (optional)
%   model       (optional)
%   schedule    (optional)
%   initState   (optional)
%   wellSols    (optional, only if simulation executed)
%   states      (optional, only if simulation executed)
%   sim_report  (optional, only if simulation executed)


% SEE ALSO:
%   runSleipner, analyseSleipner (co2lab/experimental/project/tests/)

%%
moduleCheck('co2lab','ad-core','opm_gridprocessing','mex','deckformat', ...
    'coarsegrid','upscaling','incomp','mrst-experimental');
mrstVerbose on
gravity on;


%%
opt = struct(   'gridname',             'IEAGHGmodel', ...
                'refineLevel',          1, ...
                'extractSubgrid',       false, ...
                'addPert2Grid',         false, ...
                'pertAmp',              [], ...
                'ratecase',             'SPE', ...
                'num_years',            [], ...
                'modifyPerm',           false, ...
                'modifyPoro',           false, ...
                'modifyRhoCO2',         false, ...
                'perm_fac',             [], ...
                'poro_fac',             [], ...
                'rhoCO2_fac',           [], ...
                'testSolverOptions',    false, ...
                'useSensModel',         false, ...
                'runSimulation',        true, ...
                'askBeforeSolver',      true, ...
                'askBeforeSaving',      true, ...
                'dirName',              'SleipnerSims'       );
            
opt = merge_options(opt, varargin{:});
    

    %% Load grid with given refinement/coarsening specified (if any):
    % makeSleipnerModelGrid() looks for file or generates it from grdecl files
    % and writes .mat file. Output is the variables G, Gt, rock, rock2D.
    % Grid options: 'ORIGINALmodel', 'IEAGHGmodel', 'INHOUSEmodel'
    % Coarsening level specified using -2, -3, etc.
    fprintf(['\nYou have chosen to use the ' opt.gridname ' grid.\n'])
    fprintf(['You have chosen to refine the model grid ',num2str(opt.refineLevel),' times.\n'])
    fprintf('\nGetting grid...\n\n')
    [ G, Gt, rock, rock2D ] = makeSleipnerModelGrid('modelName', opt.gridname, 'refineLevel',opt.refineLevel);
    fprintf('\n\nGrid obtained.\n')
    
    if opt.extractSubgrid
        % Extract subgrid from grid:
        dx = (var.wellXcoord-Gt.cells.centroids(:,1));
        dy = (var.wellYcoord-Gt.cells.centroids(:,2));
        ind = dx<1e3 & dx>-1.5e3 & dy<1.3e3 & dy >-4e3;
        %clf,plotCellData(G,double(ind))
        G = removeCells(G,find(~ind));
        Gt = topSurfaceGrid(G);
        % Get perm and poro of subgrid:
        rock2D.perm = rock2D.perm(ind,:);
        rock2D.poro = rock2D.poro(ind);
    end
    
    if opt.addPert2Grid
        % Introduce a perturbation in top surface:
        % TODO: implement more options for perturbation (as function?)
        if ~isempty(opt.pertAmp)
            % Adjust z values (comes from resUtsira.m)
            Gt_adjusted = Gt;
            %adjust = @(z, amp) z + amp*(2*(rand(size(z)) - 0.5));
            %adjust = @(z, amp) z + ( -amp + (amp-(-amp)).*rand(size(z)) ); % rand number between -amp and +amp
            %adjust = @(z, amp1) z + G
            %Gt_adjusted.cells.z = adjust(Gt_adjusted.cells.z, opt.pertAmp);
            %Gt_adjusted.faces.z = adjust(Gt_adjusted.faces.z, opt.pertAmp);
            %Gt_adjusted.nodes.z = adjust(Gt_adjusted.nodes.z, opt.pertAmp);
            midd=mean(Gt.nodes.coords,1);
            dd=bsxfun(@minus,Gt_adjusted.nodes.coords,midd);
            dd2=sum(dd.^2,2);
            dd2max=max(dd2);
            Gt_adjusted.nodes.z=Gt_adjusted.nodes.z+opt.pertAmp*exp(-40*dd2/dd2max);
            % Recompute geometry to get correct centroids
            Gt_adjusted = computeGeometryVE(Gt_adjusted);
            
        else
            fprintf('\nYou must specify perturbation amplitude\n.')
        end
        % figure; plotCellData(Gt, Gt.cells.z - Gt_adjusted.cells.z);
        % view(3); colorbar
        Gt = Gt_adjusted;
    end
    
    
    %% CO2 entry into layer 9
    % Physical coordinate of injection well (Singh et al. 2010)
    var.wellXcoord      = 4.38516e5;
    var.wellYcoord      = 6.47121e6;
    
    % Index of the closest cell to the physical well location
    dv              = bsxfun(@minus, Gt.cells.centroids(:,1:2), [var.wellXcoord, var.wellYcoord]);
    [v,i]           = min(sum(dv.^2, 2));   
    wellCellIndex   = i; % or Gt.cells.indexMap(i);
    [i, j]          = ind2sub(Gt.cartDims, wellCellIndex);
    % Cartesian coordinate that wellCellIndex corresponds to:
    var.wellCoord_x = Gt.cells.centroids(wellCellIndex,1);
    var.wellCoord_y = Gt.cells.centroids(wellCellIndex,2);
    var.wellCoord_z = 0;
    
    

    
    
    %% Specify fluid properties:
    [rho, mu, sr, sw]   = getValuesSPE134891();
    water_density       = rho(1) * kilogram / meter ^3;
    rhoCref             = rho(2) * kilogram / meter ^3;
    
    seafloor_temp       = 7; % Celsius
    seafloor_depth      = 100; % meters
    temp_gradient       = 35.6; % Celsius / km
    
    caprock_temperature = 273.15 + seafloor_temp + ...
        (Gt.cells.z - seafloor_depth) / 1e3 * temp_gradient;    % Kelvin
    
    
    % initial state
    initState.pressure  = Gt.cells.z * norm(gravity) * water_density;   % hydrostatic pressure, in Pa=N/m^2
    initState.s         = repmat([1 0], Gt.cells.num, 1);               % sat of water is 1, sat of CO2 is 0
    initState.sGmax     = initState.s(:,2);                             % max sat of CO2 is initially 0
    initState.rs        = 0 * initState.sGmax;                          % initially 0
    
    
    % reference reservoir conditions:
    % i.e., pressure at which water density equals the reference density (a
    % reference for linear compressibilities)
    ref_p               = mean(initState.pressure); 
    ref_t               = mean(caprock_temperature); % Kelvin
    
    %[ water_compr_val, co2_compr_val, pvMult ] = getCompressValues();
    water_compr_val     = 4.3656e-10/barsa;
    co2_compr_val       = 1.2388e-8/barsa;
    pvMult              = 1e-5/barsa; % Odd will check on this value
    
    isDissOn            = false;
    dis_max             = (53 * kilogram / meter^3) / rhoCref; % from CO2store

    % kwm? 0.75, 0.54 in Appendix of Singh et al 2010.

    
    %% Parameter modification (optional)
    % Use default parameter modifier factors if they are unspecified:
    if (opt.modifyPerm && isempty(opt.perm_fac))
        opt.perm_fac    = 3;
    end
    if (opt.modifyPoro && isempty(opt.poro_fac))
        opt.poro_fac    = 0.6;
    end  
    if (opt.modifyRhoCO2 && isempty(opt.rhoCO2_fac))
        opt.rhoCO2_fac  = 2/3;
    end
    % Perform modification:   
    if opt.modifyPerm
        disp('Original rock permeabilities are being modified ...')
        rock.perm   = rock.perm .* opt.perm_fac;
        rock2D.perm = rock2D.perm .* opt.perm_fac;
    end
    if opt.modifyPoro
        disp('Original rock porosities are being modified ...')
        rock.poro   = rock.poro .* opt.poro_fac;
        rock2D.poro = rock2D.poro .* opt.poro_fac;
    end
    if opt.modifyRhoCO2
        disp('Original CO2 density value is being modified ...')
        rhoCref = rhoCref * opt.rhoCO2_fac; 
    end

    
    %% Create fluid:
    fluid = makeVEFluid(Gt, rock2D, 'sharp interface', ...
        'fixedT'      , caprock_temperature, ...
        'wat_mu_ref'  , mu(1), ...
        'co2_mu_ref'  , mu(2), ...
        'wat_rho_ref' , water_density, ...
        'co2_rho_ref' , rhoCref, ...
        'wat_rho_pvt' , [water_compr_val, ref_p], ...
        'co2_rho_pvt' , [co2_compr_val, ref_p], ...
        'pvMult_p_ref', ref_p, ...
        'pvMult_fac'  , pvMult, ...
        'residual'    , [sw, sr] , ...
        'dissolution' , isDissOn, 'dis_max', dis_max);
    

    %% Schedule: set up time steps, well controls, and get discrete rates
   
    % First get time step sizes and associated control well numbers:
    [ schedule, var_tmp ] = ...
        getInjRatesAndTimeSteps( 'ratecase',opt.ratecase, 'num_years',opt.num_years );
    % schedule only contains step.val and step.control so far
    
    % add var_tmp to current var
    var.inj_year    = var_tmp.inj_year;
    var.inj_rates   = var_tmp.inj_rates;
    % var now contains inj_rates and inj_year


    % Then add control wells:
    for i = 1:numel(var.inj_rates)
        schedule.control(i).W = addWell([], Gt.parent, rock2D, ...
            wellCellIndex, 'name', sprintf('W%i', i), 'Type', 'rate', ...
            'Val', var.inj_rates(i), 'comp_i', [0 1]);
            % inj_rate should be mass rate / fluid.rhoGS
    end


    %% Schedule: add boundaries
    
    % First get the faces of the boundaries. face.neighbors are the indices
    % of the cells on either side of the faces, i.e., face.neighbor(100,1)
    % and face.neighbor(100,2) give the index of the cells on either side
    % of face with index 100. Any 0 cell index means there is no cell,
    % i.e., the face is along an external boundary of the domain. Thus
    % bdryFaces may be obtained by finding all the face indices that
    % contain a 0 cell index on either side.
    bdryFaces   = find( Gt.faces.neighbors(:,1).*Gt.faces.neighbors(:,2) == 0 );
    % same as bdryFaces = any(Gt.faces.neighbors==0,2);
    bdryVal     = Gt.faces.z(bdryFaces) * water_density * norm(gravity);
    bc          = addBC( [], bdryFaces, 'pressure', bdryVal, 'sat', [1 0] );   
    for i = 1:numel(schedule.control)
        schedule.control(i).bc = bc;
    end


    %% Solver testing:
    if opt.testSolverOptions
        % TODO: implement following options as varargin, or simply test
        % solver options by setting parameters here:
        %elipticsolver = BackslashSolverAD();
        %elipticsolver = AGMGSolverAD('tolerance',1e-3,'reuseSetup',false);
        %linearsolver = CPRSolverAD('tolerance',1e-5,'ellipticSolver',elipticsolver);
        linearsolver = BackslashSolverAD();
        it_target=5;
        rampup_dt = 20*day;
        dtmax=1*year;
        timestepper = ...
           IterationCountTimeStepSelector('targetIterationCount',  it_target,...
                                          'minRelativeAdjustment', sqrt(eps),...
                                          'maxRelativeAdjustment', 4,...
                                          'firstRampupStep',       rampup_dt,...
                                          'verbose', true,'maxTimestep',dtmax);
        
        nonlinearsolver = NonLinearSolver('timeStepSelector', timestepper, ...
                                'maxiterations', 2*it_target,'LinearSolver',linearsolver);
    
        nonlinearsolver = NonLinearSolver('LinearSolver',linearsolver);
        model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);
        tmpschedule=schedule;
        tmpschedule.step.control=tmpschedule.step.control(1);
        tmpschedule.step.val=tmpschedule.step.val(1);

        profile off
        profile on
        [wellSols, states, sim_report] = simulateScheduleAD(initState, model, tmpschedule,'NonLinearSolver',nonlinearsolver);
        profile off;
        profile viewer
        figure(),plotCellData(Gt,states{end}.s(:,1))
    end

    %% Create model:
    if opt.useSensModel
        % make sensitivity model
        clear CO2VEBlackOilTypeModelSens
        smodel = CO2VEBlackOilTypeModelSens(Gt, rock2D, fluid);
        % add factors to initState and schedule.control
        initState.dz      = zeros(Gt.cells.num,1);
        initState.rhofac  = 1;
        initState.permfac = 1;
        initState.porofac = 1;
        for i = 1:numel(schedule.control)
            schedule.control(i).dz      = zeros(Gt.cells.num,1);
            schedule.control(i).rhofac  = 1;
            schedule.control(i).permfac = 1;
            schedule.control(i).porofac = 1;
        end
        % rename variable
        model = smodel;
%         % run simulation using sensitivity model
%         [wellSols, states, sim_report] = simulateScheduleAD(initState, smodel, schedule);
%         states = addHeightData(states,Gt,fluid);
    else
        model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);
    end

    %% Run simulation:
    % Option to allow for function to run but break before the simulation.
    % This allows for the model set-up and inputs to be analyzed without
    % having to run a simulation.
    if ~opt.runSimulation
        disp('Simulation was not run.')
    else
        % If ask=true, simulation is executed after obtaining user consent.
        % If ask=false, no user consent is obtained and simulation is run
        % automatically.
        if (opt.askBeforeSolver && userConsent('Do you want to proceed to solver?'))
            [wellSols, states, sim_report] = simulateScheduleAD(initState, model, schedule);
            saveResults( opt, var, model, schedule, wellSols, states, sim_report)

        elseif ~opt.askBeforeSolver
            [wellSols, states, sim_report] = simulateScheduleAD(initState, model, schedule);
            saveResults( opt, var, model, schedule, wellSols, states, sim_report)

        end
    end

    %% Variables to pass out of function (optional):
    if nargout ~= 0
        varargout{1} = opt;
        varargout{2} = var;
        varargout{3} = model;
        varargout{4} = schedule;
        varargout{5} = initState;
        if nargout > 5
            if exist('wellSols','var')
                varargout{6} = wellSols;
            else
                varargout{6} = [];
            end
            if nargout > 6
                if exist('states','var')
                    varargout{7} = states;
                else
                    varargout{7} = [];
                end
                if nargout > 7
                    if exist('sim_report','var')
                        varargout{8} = sim_report;
                    else
                        varargout{8} = [];
                    end
                end
            end
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%       HELPER FUNCTIONS      %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [ cw, cg, cr ] = getCompressValues()
        
        % If .mat files are not available, the values which are computed
        % are as follows:
        % water_compr_val     = 4.3656e-10/barsa;
        % co2_compr_val       = 1.2388e-8/barsa;
        % pvMult              = 1e-5/barsa; % Odd will check on this value
        
        a_w = SampledProp2D('rho','Water_100000_400000000_278_524_800_800_D.mat');
        a_g = SampledProp2D('rho','CarbonDioxide_100000_400000000_278_524_800_800_D.mat');
        
        cw = a_w.rhoDP(ref_p,ref_t) / a_w.rho(ref_p,ref_t);
        cg = a_g.rhoDP(ref_p,ref_t) / a_g.rho(ref_p,ref_t);
        cr = [];
        
    end

    % ---------------------------------------------------------------------




    %%%%%%%%%%%%%%%%%%%%%   END OF HELPER FUNCTIONS    %%%%%%%%%%%%%%%%%%%%

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        OTHER FUNCTIONS       %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveResults( opt, var, model, schedule, wellSols, states, sim_report )

    % Note: only variables passed into this function will be saved to .mat
    % file

    % first, close all figures
    close all
    %if strcmpi(opt.mycase(1:3), 'use') && strcmpi(opt.mycase(end-5:end), '_model')
        %name = opt.gridname(1:end-5);
    %end
    %if strcmpi(opt.ratecase(1:3), 'use')
        %rateName = opt.ratecase;
    %end
    [modstr1, modstr2, modstr3] = deal('_');
    if opt.modifyPerm
        modstr1 = ['_PermFac', num2str(opt.perm_fac), '_'];
    end
    if opt.modifyPoro
        modstr2 = ['PoroFac', num2str(opt.poro_fac), '_'];
    end
    if opt.modifyRhoCO2
        modstr3 = ['RhoFac', num2str(opt.rhoCO2_fac), '_'];
    end
    modstrs = [modstr1 modstr2 modstr3];
    fileName = [opt.gridname(1:end-5) 'refNum' num2str(opt.refineLevel) '_' opt.ratecase 'rates' modstrs datestr(clock,30) '.mat'];

    % clear some unneeded variables before saving
    clearvars modstrs modstr1 modstr2 modstr3

    % To avoid automatically saving simulation results which can create large
    % .mat files, user consent is required if ask=true. Otherwise, if
    % ask=false, saving occurs automatically.
    if (opt.askBeforeSaving && userConsent('Do you want to save simulation results?'))
        mkdir(opt.dirName)
        save([opt.dirName '/' fileName],'-v7.3');

    elseif ~opt.askBeforeSaving
        mkdir(opt.dirName)
        save([opt.dirName '/' fileName],'-v7.3');
    end
end

