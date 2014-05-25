function runVEScenario(comp_model, scenariofile, savename, square_domain)
%
%
% SYNOPSIS:
%   function runVEScenario(model, savename)
%
% DESCRIPTION:
%
% PARAMETERS:
%   comp_model - 1 - incompressible; 2 - semi-compressible; 3 - fully compressible
%   savename   - file in which to save the results
%
% RETURNS:
%   runVEScenario(model, - 
%
% EXAMPLE:
%
% SEE ALSO:
%
    moduleCheck('ad-fi', 'ad-refactor');
    gravity on;
    if ~exist('square_domain')
        % Default is a single row of cells
        square_domain = false;
    end
    

    %% load scenario-specific parameters
    refcell_ix = 1; % reference cell for computing constant density - can be
                    % changed in scenario file 
    
    run(scenariofile);
    
    %% define boundary conditions
    pCap = hydrostaticPressure(Gt, rhoW, 1*atm, slope, slopedir, zeros(Gt.cells.num, 1));
    bc   = pside([], Gt, 'LEFT',  ...
                 pCap(sum(Gt.faces.neighbors(boundaryFaceIndices(Gt,'LEFT'), :), 2)));  
    bc   = pside(bc, Gt, 'RIGHT', ...
                 pCap(sum(Gt.faces.neighbors(boundaryFaceIndices(Gt, 'RIGHT'),:),2)));
    if square_domain
        % assert(slope==0);
        % bc = pside(bc, Gt, 'FRONT', ...
        %            pCap(sum(Gt.faces.neighbors(boundaryFaceIndices(Gt,'FRONT'), :), 2)));  
        % bc = pside(bc, Gt, 'BACK' , ...
        %            pCap(sum(Gt.faces.neighbors(boundaryFaceIndices(Gt,'BACK'), :), 2)));  
    end
    
    %% Define model
    CO2obj = CO2props('rho_big_trunc','');
    depth = computeRealDepth(Gt, slope, slopedir, 0); % slope-adjusted depth

    switch comp_model 
      case 1
        % determine constant rho value
        tref = ref_temp + depth(refcell_ix) * temp_grad/1000;
        rhoC = CO2obj.rho(pCap(refcell_ix), tref);
        % make the EOS return only this value of rho, regardless of P and T
        CO2obj.rho = @(p, t) rhoC * ones(numel(double(p)), 1);
        CO2obj.compressible = 'incompressible';
      case 2
        CO2obj.compressible = 'horizontal';
      case 3
        CO2obj.compressible = 'full';
        CO2obj = fullCompressibleCO2BrineModel.defineAdditionalDerivatives(CO2obj);
      otherwise 
        error('Incorrect choice of compressibility model provided.');
    end
    model = fullCompressibleCO2BrineModel(Gt, rock, tinfo, ...
            'EOSCO2', CO2obj, ...
            'constantVerticalDensity', strcmp(CO2obj.compressible, 'horizontal'), ...
            'rhoBrine', rhoW, ...
            'mu', mu, ...
            'slope', slope, ...
            'slopedir', slopedir);

    %% Initialize state
    % @@ NB: computation of hydrostatic pressure below must be changed if
    % reservoir is tilted
    
    mass_ixs = h0<0;
    approximate = true;
    h0(mass_ixs) = convertMassesToHeight(-h0(mass_ixs), CO2obj, ...
                                         pCap(mass_ixs), tinfo, rhoW, slope, ...
                                         depth(mass_ixs), ...
                                         Gt.cells.volumes(mass_ixs), ...
                                         rock.poro(mass_ixs), approximate);
    pI = hydrostaticPressure(Gt, rhoW, 1*atm, slope, slopedir, h0);
    state = struct('pressure', pI, 'h', h0, 'wellSol', []);
    if ~isempty(schedule.W)
        state.wellSol = struct('bhp', pI(1), ...
                               'qGs', schedule.W(1).val);
    end
    
    %% Run schedule
    [wellSols, states] = runScheduleRefactor(state, model, schedule, 'bc', bc);
    
    %% compute caprock values and save result
    quick = true;
    for i = 1:numel(states)
        states{i} = model.includeComputedCaprockValues(states{i}, quick);
    end        

    result = makeResultStructure(states, Gt, rock);
    save(savename, 'result', 'Gt');
end



