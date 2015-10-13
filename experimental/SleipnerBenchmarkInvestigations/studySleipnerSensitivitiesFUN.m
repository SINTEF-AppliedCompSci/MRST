function [ varargout ] = studySleipnerSensitivitiesFUN( initState, smodel, schedule, wellSols, states, varargin )
% Study sensitivties and optimization of properties (such as top surface
% elevation, CO2 density, CO2 entry rates, etc) related to simulation
% matching of the Sleipner Benchmark.

% Assess sensitivity of Inputs to the match between simulated and observed CO2 heights.
    % Inputs include model grid, CO2 density, porosity, permeability.


% SYNOPSIS:
%
% studySleipnerSensitivitiesFUN( model, states )
% studySleipnerSensitivitiesFUN( model, states, 'plumes_base', states_base )



% PARAMETERS:
%   model       - the model used in the simulation scenario which produced states.
%   states      - the so called 'simulated' plumes.
%   plumes_base - the so called 'observed' plumes. If not passed into
%                 function, plumes_base is loaded from observed plume .mat
%                 files.


% OPTIONS:


% RETURNS:
%   dobj_dz - total derivative of objective function (L2-norm) wrt height of top surface



% SEE ALSO:
%   studySleipnerBenchmarkFUN.m


%%

    opt.plumes_base = [];
    
    opt = merge_options(opt, varargin{:});
    
    if isempty(opt.plumes_base)
        % default is to use the observed CO2 plume outlines. The outlines
        % are loaded, and then the CO2 heights inside the plume outlines
        % are computed, wrt a top-surface grid. Top-surface grid is taken
        % to be the same as the one passed in that corresponds to states.
        fprintf('\n Using observed CO2 plume outlines as plumes_base. \n')
        opt.plumes_base = getLayer9CO2plumeOutlines();
        
        fprintf('\n Computing observed CO2 heights wrt top surface of smodel.G. \n')
        [opt.plumes_base, ~, ~, ~] = makeSurfaceDataAndPlots(opt.plumes_base, smodel.G);
    else
        % ensure plumes_base.h exists
        fprintf('\n Using your input argument as plumes_base. \n')
        for i = 1:numel(opt.plumes_base)
            assert(isfield(opt.plumes_base{i},'h'), ...
            'CO2 heights inside plume outline have not been computed yet.')
        end
        
    end


    moduleCheck('co2lab','ad-core','opm_gridprocessing','mex','deckformat', ...
    'coarsegrid','upscaling','incomp','mrst-experimental','optimization');
    mrstVerbose on
    gravity on
    
    % For generating plots
    %myplotCellData = @(G, data) plotCellData(G, data, 'EdgeColor','none');
    
    
    %% Renaming variables
    Gt          = smodel.G;
    fluid       = smodel.fluid;
    newplumes   = opt.plumes_base;
    

    
    %% Grid Sensitivity:
    % Using functions topsurface and topfit, get elevations of top
    % corresponding to model grid and a planar surface (which is fitted to
    % plume outline). Get difference between these two 'tops' = hCO2. Get
    % difference between 'observed' (hCO2) and 'simulated' CO2 height.
    
    %%[ plumes, topsurface, topfit, hCO2 ] = makeSurfaceData(plumes, Gt);
    %[ plumes, topsurface, topfit, hCO2 ] = makeSurfaceDataAndPlots(newplumes, Gt);
	%[ plumes_comb ] = plotDiff( Gt, plumes, sim_report, states, topsurface );
     

    %a = matchToData(model, wellSols, states, schedule, newplumes);

    
    %% Compute gradients using adjoint method:
    
    % Objective function: the match is quantified by calculating the
    % L2-norm of 'observed' to 'simulated' CO2 heights.
    obj_funs = @(wellSols, states, schedule, varargin) matchToDataSens(smodel, wellSols, states, schedule, newplumes, varargin{:});
    objhs = @(tstep) obj_funs(wellSols, states, schedule, 'ComputePartials', true, 'tStep', tstep);
    
    % Gradients
    gsc   = computeGradientAdjointAD(initState, states, smodel, schedule, objhs, 'ControlVariables', {'well','scell','mult'});

    % Get total derivatives
    % (i.e., sum each partial derivatives over all time steps)
    dobj    = cell(size(gsc,1),1);
    dobj{1} = zeros(1,1);               % partial derivatives wrt well rates
    dobj{2} = zeros(Gt.cells.num,1);    % partial derivatives wrt height of top surface
    dobj{3} = zeros(3,1);               % partial derivatives wrt global multiplies: CO2 density, perm, poro
    for i = 1:size(gsc,2)
        for j = 1:size(gsc,1);
            dobj{j} = dobj{j} + gsc{j,i};
        end
    end
    dobj_dz = dobj{2}; % total derivative wrt height of top surface

    
    %% Variables to pass out of function
    if nargout ~= 0
        varargout{1} = dobj_dz;
    end
    
     
     
end

