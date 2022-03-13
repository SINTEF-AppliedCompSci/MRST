function [updata, report] = upRelPerm(block, updata, method, varargin)
% Upscaling of relative permeability

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opt = struct(...
    'nsat',        20,         ... % Number of upscaled sat. values
    'values',      [],         ... % Specify the values
    'dims',        1:3,        ... % What dimensions to upscale
    'dp',          1*barsa,    ... % Pressure drop
    'savesat',     false,      ... % Save saturation distributions
    'absmethod',   'pressure', ... % One-phase upscaling method
    'verbose',     false       ... % Print progress to console
    );
opt = merge_options(opt, varargin{:});

wantReport  = nargout > 1;
wantSatDist = opt.savesat && wantReport;
timeStart   = tic;

dims  = opt.dims;
ndims = length(dims);
if ~isempty(opt.values), opt.nsat = numel(opt.values); end
nvals = opt.nsat;
rock  = block.rock;
fluid = block.fluid;

assert(isa(block, 'GridBlock'));
assert(isfield(updata, 'perm'), 'One phase upscaling must be run first.');
assert(~isempty(method), 'Method must be set');

% Get upscaled absolute permeability
Kup = updata.perm;
if numel(Kup)>ndims && numel(Kup)==3
    % We are given upscaled absperm in all three dimensions
    % Extract the dimensions we need
    Kup = Kup(dims);
end
assert(numel(Kup)==ndims);

% Pore volume
pvTot = sum(block.pv);

% Allocate space
krW = cell(1, ndims);
krO = cell(1, ndims);
if wantSatDist
    satdist = cell(nvals,1);
    if strcmpi(method, 'capillary-viscous-dist')
        satdistff = cell(nvals,1);
    end
end

if ~isempty(opt.values)
    % The values are explicity specified
    values = opt.values;
else
    % Get input values depending on the method
    values = getValues(block, updata, method, nvals);
end

dispif(opt.verbose, '\nStarting relperm upscaling\n');
dispif(opt.verbose, '   #      sW   Time(s)\n');
start = tic;

% Loop over the input values
for iv = 1:nvals
    
    startSatValue = tic;
    
    %dispif(opt.verbose, 'Saturation value %d of %d.\n', iv, nvals);
    
    val = values(iv);
    
    % Get saturation distribution for this input value
    [sW0, sWff] = valueDistribution(block, method, val, opt.savesat);
    
    % Save if requested
    if wantSatDist
        satdist{iv} = sW0;
        if strcmpi(method, 'capillary-viscous-dist')
            satdistff{iv} = sWff;
        end
    end
    
    % Loop over the dimension
    for id = 1:ndims
        d = dims(id); % Current dimension
        
        % Update the state for this dimension
        if isempty(opt.values) && ( iv==1 || iv==nvals )
            % If the values are not specified from outside, we know that at
            % the end points, nothing changes, and we do not need to runa
            % the simulation to steady state.
            sW = sW0;
        else
            sW = directionDistribution(block, method, sW0, d, opt.dp, ...
                opt.verbose);
        end
        
        % Upscaled saturation
        sWup = sum(sW.*block.pv)./pvTot;
        
        % Relative phase permeability
        if isfield(fluid, 'krOW')
            kr{1} = fluid.krOW(1-sW);
        else
            kr{1} = fluid.krO(1-sW);
        end
        kr{2} = fluid.krW(sW);
        
        % Set permeability as K*kr
        rock_Kkr = rock;
        
        
        % Loop over phases
        for p = 1:2
            
            % If the permeability is anisotropic, we multiply each
            % direction with the relperm kr
            rock_Kkr.perm = bsxfun(@times, rock.perm, kr{p});
            
            if all(all(rock_Kkr.perm == 0))
                % We will get no fluid motion
                krKU = 0;
            else
                if any(any(rock_Kkr.perm < 0))
                    error('Some pseudo perm values are negative!');
                end
                if any(any(rock_Kkr.perm == 0))
                    % To avoid a singular matrix when performing one
                    % phase upscaling, the zero permeabilities are set
                    % to something larger than zero.
                    ep = eps(mean(mean(...
                        rock_Kkr.perm(rock_Kkr.perm>0)) ));
                    rock_Kkr.perm(rock_Kkr.perm == 0) = ep*1e1;
                end
                
                % Perform one phase upscaling with the altered
                % permeability field
                block.rock = rock_Kkr;
                data = upAbsPerm(block, [], 'dims', d, ...
                    'dp', opt.dp, 'method', opt.absmethod);
                krKU = data.perm;
            end
            
            % Compute upscaled relperm value
            krup = krKU / Kup(id);
            
            if p==1
                krO{id}(iv,:) = [sWup, krup];
            else
                krW{id}(iv,:) = [sWup, krup];
            end
            
        end
        
        
    end % End of dimension loop
    
    secs = toc(startSatValue);
    dispif(opt.verbose, ' %2d/%d   %1.2f   %1.2f\n', ...
        iv, nvals, sWup, secs);
    
end % End of input value loop

secs = toc(start);
dispif(opt.verbose, ' Total time:   %4.2f s\n', secs);

% Check for upscaled values outside range. We simply force the values
% inside valid range.
outside = false;
for id = 1:ndims
    inx=krO{id}(:,2)>1; krO{id}(inx,2)=1; if any(inx),outside=true; end
    inx=krO{id}(:,2)<0; krO{id}(inx,2)=0; if any(inx),outside=true; end
    inx=krW{id}(:,2)>1; krW{id}(inx,2)=1; if any(inx),outside=true; end
    inx=krW{id}(:,2)<0; krW{id}(inx,2)=0; if any(inx),outside=true; end
end

% Store upscaled data to structure
updata.krO = krO;
updata.krW = krW;

totalTime = toc(timeStart);
if wantReport
    report.method  = method;
    report.dims    = dims;
    report.nsat    = nvals;
    report.dp      = opt.dp;
    report.time    = totalTime;
    report.valsOutsideRange = outside;
    if wantSatDist
        report.satdist = satdist;
        if strcmpi(method, 'capillary-viscous-dist')
            report.satdistff = satdistff;
        end
    end
end


end


%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------

function values = getValues(block, updata, method, nvals)

switch method
    case 'flow'
        values = linspace(0, 1, nvals)';
        
    case {'capillary', 'capillary-viscous-dist'}
        assert(isfield(updata, 'pcOW'), ...
            'Run capillary curve upscaling first');
        swUMin = updata.pcOW(1,1);
        swUMax = updata.pcOW(end,1);
        sW     = linspace(swUMin, swUMax, nvals)';
        values = interp1(updata.pcOW(:,1), updata.pcOW(:,2), sW);
        
    case 'viscous'
        ffdata = upFracFlowOW(block, []);
        swUMin = ffdata.ffW(1,1);
        swUMax = ffdata.ffW(end,1);
        sW     = linspace(swUMin, swUMax, nvals)';
        values = interp1(ffdata.ffW(:,1), ffdata.ffW(:,2), sW);
        
    case {'capillary_grav', 'capillary-viscous-dist_grav'}
        assert(isfield(updata, 'pcOW_bot'), ...
            'Run gravity capillary curve upscaling first');
        swUMin = updata.pcOW_bot(1,1);
        swUMax = updata.pcOW_bot(end,1);
        sW     = linspace(swUMin, swUMax, nvals)';
        values = interp1(updata.pcOW_bot(:,1), updata.pcOW_bot(:,2), sW);
        
    otherwise
        error(['Method ''' method ''' not recognized.']);
end

end


function [sW, sWff] = valueDistribution(block, method, val, savesat)
% Get the saturation distribution for the current value, depending on the
% method chosen. The returned saturation may be updatad depending on the
% direction later.

sWff  = [];
G     = block.G;
fluid = block.fluid;

% Get saturation values
switch method
    case 'flow' % Flow based simulation
        assert(all(isfield(fluid, {'swir','sor','satnum'})), ...
            'Fluid structure does not have necessary fields.');
        
        sW = val;
        
        % We have different saturation regions. The sW values are
        % mapped from [0 1] to the different [swir 1-sor] intervals.
        swir = fluid.swir;
        sor  = fluid.sor;
        
        % sW mapped from scalar to vector of length nRegions
        sW = swir + sW*( (1-sor) - swir );
        
        % Map sW values to each cell
        sW = sW(fluid.satnum);
        
    case 'capillary' % Capillary limit
        assert(isfield(fluid, 'pcOWInv'), ['To use capillary limit '...
            'upscaling, the fluid must have field ''pcOWInv''.']);
        pcOW = val;
        sW   = fluid.pcOWInv(pcOW.*ones(G.cells.num,1));
        
    case 'viscous' % Viscous limit
        assert(isfield(fluid, 'fracFlowInv'), ['To use viscous limit '...
            'upscaling, the fluid must have field ''fracFlowInv''.']);
        ff   = val;
        sW   = fluid.fracFlowInv(ff.*ones(G.cells.num,1));
        
    case 'capillary_grav' % Capillary limit with gravity
        % The pcOW values are given at some chosen datum
        assert(isfield(fluid, 'pcOWInv'), ['To use capillary limit '...
            'upscaling, the fluid must have field ''pcOWInv''.']);
        pcOW = val; % pcOW at some level
        
        if isfield(fluid, 'rhoO')
            rhoO = fluid.rhoO;
        else
            rhoO = fluid.rhoOS;
        end
        if isfield(fluid, 'rhoW')
            rhoW = fluid.rhoW;
        else
            rhoW = fluid.rhoWS;
        end
        dRho = rhoW - rhoO;
        g    = 9.8066; % HARDCODED
        zi   = max(G.cells.centroids(:,3)) - G.cells.centroids(:,3);
        grav = -dRho.*g.*zi; % NOTE: Not sure about sign here
        
        sW   = fluid.pcOWInv(pcOW - grav);
        
    case 'capillary-viscous-dist' % Experimental
        
        pcOW = val;
        if savesat
            [sW, sWff] = getCapVisDist(block, pcOW);
        else
            sW = getCapVisDist(block, pcOW);
        end
        
    case 'capillary-viscous-dist_grav' % Experimental
        
        pcOW = val;
        if savesat
            [sW, sWff] = getCapVisDist(block, pcOW, 'gravity', true);
        else
            sW = getCapVisDist(block, pcOW, 'gravity', true);
        end
        
    otherwise
        error(['Method ''' method ''' not recognized.']);
end

end


function sW = directionDistribution(block, method, sW, d, dp, verbose)
% Update saturation for the current direction

% Get saturation values
switch method
    case 'flow' % Flow based simulation
        
        G = block.G;
        isPeriodic = block.periodic;
        
        % Set pressure drop
        if isPeriodic
            assert( isprop(block, 'bcp') );
            bcp = block.bcp;
            bcp.value(:) = 0;
            bcp.value(bcp.tags == d) = dp;
            bc  = [];
        else
            % Cells neighbouring the outer faces
            nc = @(i) sum(block.G.faces.neighbors(block.faces{d}{i},:),2);
            
            sWf = sW(nc(1));
            bc = addBC([], block.faces{d}{1}, ...
                'pressure', dp, 'sat', [sWf 1-sWf]);
            
            sWf = sW(nc(2));
            bc = addBC(bc, block.faces{d}{2}, ...
                'pressure',  0, 'sat', [sWf 1-sWf]);
            
            bcp = [];
        end
        
        % Solve steady state flow
        state1 = simulateToSteadyStateADI(G, block.rock, ...
            block.fluid, sW, 'bc', bc, 'bcp', bcp);
        
        % Extract saturation
        sW = state1.s(:,1);
        
    case 'capillary' % Capillary limit
        % Do nothing.
    case 'viscous' % Viscous limit
        % Do nothing.
    case 'capillary_grav' % Capillary limit with gravity
        % Do nothing.
    case 'capillary-viscous-dist' % Experimental
        % Do nothing.
    case 'capillary-viscous-dist_grav' % Experimental
        % Do nothing.
    otherwise
        error(['Method ''' method ''' not recognized.']);
end

end





