function [ schedule ] = setSchedule_extras( Gt, rock2D, wcells, wtype, ...
                                isteps, itime, msteps, mtime, ...
                                varargin)
% Put well controls and time steps into schedule. Wells can be either
% pressure-controlled ('bhp') or rate-controlled ('rate'). Wells in
% migration period are assigned as rate-controlled, with rates = minval.

% @@ implement to handle production wells.

    opt.initOverP   = 1 * mega * Pascal; % used if wtype = 'bhp'
    opt.wellradius  = 0.3;
    opt.minval      = 0; % sqrt(eps)    % minimum rate for migration period
    opt.wqtots      = [];               % required if wtype = 'rate'
    opt.initState   = [];               % required if wtype = 'bhp'
    %opt.single_control = [];
    
    opt = merge_options(opt, varargin{:});
   
    assert(isteps>0);
    if msteps == 1
        msteps = 2;
        warning(['If migration is happening, we need at least two steps, in ' ...
                 'order to have a transition step.  Migration step has been ' ...
                 'increased to two.']);
    end

    
    %allWells = [opt.wellCoords_inj; opt.wellCoords_prod];
    
    numWells = size(wcells,1);
    
    if strcmpi(wtype,'rate')
        
        assert(~isempty(opt.wqtots));
        assert(numel(opt.wqtots) == numWells);
        % computing fixed rates
        wrates = opt.wqtots / itime;
        wvals = wrates; % %[ opt.inj_rate; -opt.prod_rate ];
        wvals(wvals==0) = opt.minval;
        
        
    elseif strcmpi(wtype,'bhp')

        assert(~isempty(opt.initState));
        % bhp is computed using initial pressure of formation plus the
        % maximum possible overpressure specified by user. Instead of
        % limiting the bhp by a maximum overpressure value, the overburden
        % pressure could be computed as a function of caprock depth, and
        % then the overpressure could be set to 0.9 times this overburden
        % pressure to avoid fracturing the caprock.
        wvals = opt.initState.pressure(wcells) + opt.initOverP;
        %wbhp = BHP; %ones(numWells,1) .* opt.BHP;
        assert(numel(wvals) == numWells);
        %wellVal = injBHP;
        
    else
        error('Unknown well type was specified. Options are ''rate'' or ''bhp''.')
    end
    
    W = [];
    for i = 1:numWells
        
%         % Cell index(es) of well point(s) in top grid formation
%         var.wellCellIndex(i)  = getCellIndex(Gt, allWells(i,1), allWells(i,2));
%         var.wellCoordSim(i,1) = Gt.cells.centroids(var.wellCellIndex(i),1);
%         var.wellCoordSim(i,2) = Gt.cells.centroids(var.wellCellIndex(i),2);
        
        % Assign well name (to distinguish inj vs. prod)
        if (strcmpi(wtype,'rate') && wvals(i) > 0) ...
                || (strcmpi(wtype,'bhp') && wvals(i) > 0)
            % inj
            wellname = sprintf('Winj%i', wcells(i));
            composition = [0 1];
        else
            % prod
            wellname = sprintf('Wprd%i', wcells(i));
            composition = [1 1];
        end
        
        % NB: reference depth is supplied here when adding well. If Gt.parent
        % is used to set up well, supply reference depth of
        % Gt.parents.cells.centroid(wellCellIndex,3). If Gt is used, ref depth
        % is Gt.cells.z(wellCellIndex), but it might not be necessary to supply
        % it explicitly.
        
        % Put injection/producer wells into schedule.control(1):
        % if 3D exists, addWell + convertwellsVE() is used, otherwise
        % 2D is used to create W, and then W.dZ is manually set to 0 and
        % W.WI is corrected (st it corresponds to well cell column area,
        % not area of top grid well cell)
        % NB: whether 'refDepth' is zero or otherwise, results appear to be
        % the same.
        if isfield(Gt,'parent')
            W = addWell(W, Gt.parent, rock2D, wcells(i), ...
                'name',     wellname,  ...
                'Type',     wtype, ...
                'Val',      wvals(i), ... %allRates(i), ...
                'comp_i',   composition, ...
                'Radius',   opt.wellradius, ...
                'refDepth', Gt.parent.cells.centroids(wcells(i),3) );
            %W = convertwellsVE(W, Gt.parent, Gt, rock2D, 'ip_tpf'); % do
            %outside for loop since additional field is added to W, which
            %causes number of fields in structure W not to match.
            
            %alternative: convertwellsVE_s()
            %  (however, may compute WI differently, won't contain h, and
            %  won't use the user specified 'refDepth')
            %W = convertwellsVE_s(W, Gt.parent, Gt, rock2D, 'ip_tpf');        
        else
            W = addWell(W, Gt, rock2D, wcells(i), ...
                'name',     wellname,  ...
                'Type',     wtype, ...
                'Val',      wvals(i), ... %allRates(i), ...
                'comp_i',   composition, ...
                'Radius',   opt.wellradius, ...
                'refDepth', Gt.cells.z(wcells(i)) + Gt.cells.H(wcells(i))./2 );
            W(i).dZ = 0;
            W(i).WI = W(i).WI .* Gt.cells.H(wcells(i));
        end
    end
    
    if isfield(Gt,'parent')
        W = convertwellsVE(W, Gt.parent, Gt, rock2D, 'ip_tpf');
    end
    
    
    
    % Specify well rates = 0 for control(2)
    W_shut = W;
    for i = 1:numel(W_shut)
        W_shut(i).type  = 'rate';
        W_shut(i).val   = opt.minval; 
        %W_shut(i).WI    = 0; %@@ when set to 0, solver breaks
    end

    schedule.control(1).W = W;
    schedule.control(2).W = W_shut;
    
    
    % Set up time step:
    dTi         = itime / isteps;
    dTm         = mtime / msteps;
    istepvec    = ones(isteps, 1) * dTi;
    mstepvec    = ones(msteps, 1) * dTm;
    schedule.step.val       = [istepvec; mstepvec];
    schedule.step.control   = [ones(isteps, 1); ones(msteps, 1) * 2];
    
    
%     % Constructing schedule
%     cpos = 1;
%     ctrls = [];
%     if single_control && isteps > 0 
%         schedule.control(cpos).W = W;
%         cpos = cpos+1;
%         ctrls = ones(isteps,1);
%     else
%         for i = 1:isteps
%             schedule.control(cpos).W = W;
%             cpos = cpos+1;
%         end
%         ctrls = [1:isteps]';
%     end    
%     
%     if msteps > 0
%         schedule.control(cpos).W = W;
%         for i = 1:numel(schedule.control(cpos).W)
%             schedule.control(cpos).W(i).val = opt.minval;
%         end
%         cpos = cpos+1;
%     end
% 
%     mig_ctrl = numel(schedule.control); % control step for migration, only
%                                         % relevant (and correct) if msteps > 0
% 
%     dTi = itime/isteps;
%     dTm = mtime/(msteps-1);
%     
%     istepvec = ones(isteps, 1) * dTi; 
%     mstepvec = [dTi; ones(msteps-1,1) * dTm];
%     mstepvec(2) = mstepvec(2) - dTi;
% 
%     schedule.step.val = [istepvec; mstepvec];
%     schedule.step.control = [ctrls; mig_ctrl * ones(msteps,1)];
    
    
    validateSchedule(schedule)
end

function validateSchedule(schedule)
% copied from simulateScheduleAD()

    assert (all(isfield(schedule, {'control', 'step'})));

    steps = schedule.step;

    assert (all(isfield(steps, {'val', 'control'})));

    assert(numel(steps.val) == numel(steps.control));
    assert(numel(schedule.control) >= max(schedule.step.control))
    assert(min(schedule.step.control) > 0);
    assert(all(schedule.step.val > 0));
end
