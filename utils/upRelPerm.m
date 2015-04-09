function updata = upRelPerm(block, updata, method, varargin)

opt = struct(...
    'nvalues',   20, ...
    'dims',      1:3, ...  % Dimensions to upscale
    'dp',        1*barsa ...  % Pressure drop
    );
opt = merge_options(opt, varargin{:});

dims  = opt.dims;
ndims = length(dims);
nvals = opt.nvalues;
rock  = block.rock;
fluid = block.fluid;

assert(isa(block, 'GridBlock'));
assert(isfield(updata, 'K'), 'One phase upscaling must be run first.');
assert(~isempty(method), 'Method must be set');

% Get upscaled absolute permeability
Kup = updata.K;
assert(numel(Kup)==ndims);

% Pore volume
pvTot = sum(block.pv);

% Allocate space
krW = cell(1, ndims);
krO = cell(1, ndims);

% Select columns depending on isotropic perm or not
permCol = [1 1 1];
if size(rock.perm, 2) == 3
    permCol = [1 2 3];
end

% Get input values depending on the method
values = getValues(block, method, nvals);


% Loop over the input values
for iv = 1:nvals
    
    val = values(iv);
    
    % Get saturation distribution for this input value
    sW0 = valueDistribution(block, method, val);
    
    % Loop over the dimension
    for id = 1:ndims
        d = dims(id); % Current dimension
        
        % Update the state for this dimension
        sW = directionDistribution(block, method, sW0, opt.dp);
        
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
            rock_Kkr.perm = rock.perm(:, permCol(d)) .* kr{p};

            if all(rock_Kkr.perm == 0)
                % We will get no fluid motion
                krKU = 0;
            else
                if any(rock_Kkr.perm < 0)
                    error('Some pseudo perm values are negative!');
                end

                % Perform one phase upscaling with the altered permeability
                % field
                block.rock = rock_Kkr;
                krKU = upAbsPerm(block, 'dims', d, 'dp', opt.dp);
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
    
end % End of input value loop


% Store upscaled data to structure
updata.krO = krO;
updata.krW = krW;


end


%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------

function values = getValues(block, method, nvals)

switch method
    case 'flow'
        error('TODO: Need to determine swU min and max here')
        sW = linspace(swUMin, swUMax, nvals)';
        
    case 'capillary'
        [pcFun, swUMin, swUMax] = pcVsUpscaledSw(block.G, ...
            block.rock, block.fluid);
        sW = linspace(swUMin, swUMax, nvals)';
        values = pcFun(sW); % sW upscaled -> pcOW
        
    case 'viscous'
        [ffFun, swUMin, swUMax] = fracFlowVsUpscaledSw(block.G, ...
            block.rock, block.fluid);
        sW = linspace(swUMin, swUMax, nvals)';
        values = ffFun(sW); % sW upscaled -> fractional flow
        
    case 'gravcapillary'
        [pcFun, swUMin, swUMax] = pcVsUpscaledSwGravityBinary(block.G, ...
            block.rock, block.fluid);
        sW = linspace(swUMin, swUMax, nvals)';
        values = pcFun(sW); % sW upscaled -> pcOW
        
    otherwise
        error(['Method ''' method ''' not recognized.']);
end

end


function sW = valueDistribution(block, method, val)
% Get the saturation distribution for the current value, depending on the
% method chosen. The returned saturation may be updatad depending on the
% direction later.

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
        satnum = fluid.satnum;
        sW = sW(satnum);
        
    case 'capillary' % Capillary limit
        assert(isfield(fluid, 'pcOWInv'), ['To use capillary limit '...
            'upscaling, the fluid must have field ''pcOWInv''.']);
        pcOW = val;
        sW   = fluid.pcOWInv(pcOW.*ones(G.cells.num,1));
        
    case 'viscous' % Vuscous limit
        assert(isfield(fluid, 'fracFlowInv'), ['To use viscous limit '...
            'upscaling, the fluid must have field ''fracFlowInv''.']);
        ff   = val;
        sW   = fluid.fracFlowInv(ff.*ones(G.cells.num,1));
        
    case 'gravcapillary' % Capillary limit with gravity
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
        dRho = rhoO - rhoW;
        %g    = norm(gravity);
        g    = 9.8066; % HARDCODED
        zi   = max(G.cells.centroids(:,3)) - G.cells.centroids(:,3);
        grav = dRho.*g.*zi;
        sW   = fluid.pcOWInv(pcOW - grav);
        
    otherwise
        error(['Method ''' method ''' not recognized.']);
end

end


function sW = directionDistribution(block, method, sW, dp)
% Update saturation for the current direction

% Get saturation values
switch method
    case 'flow' % Flow based simulation
        
        G = block.G;
        isPeriodic = block.periodic;
        
        % Set pressure drop
        if isPeriodic
            bcp = block.bcp;
            bcp.value(:) = 0;
            bcp.value(bcp.tags == d) = dp;
            bc  = [];
            Gp  = G;
            G   = G.parent;
        else
            bc = addBC([], block.faces{d,1}, ...
                'pressure', dp, 'sat', [sW 1-sW]);
            bc = addBC(bc, block.faces{d,2}, ...
                'pressure',  0, 'sat', [sW 1-sW]);
            Gp  = [];
            bcp = [];
        end
        
        % Initialize solution
        state0      = initResSol(G, 100*barsa, [sW, 1-sW]);
        state0.flux = zeros(G.faces.num, 2);
        
        % Solve steady state flow
        state1 = simulateToSteadyStateADI(state0, G, block.rock, ...
            block.fluid, 'bc', bc, 'Gp', Gp, 'bcp', bcp);
        
        % Extract saturation
        sW = state1.s(:,1);
        
    case 'capillary' % Capillary limit
        % Do nothing.
    case 'viscous' % Vuscous limit
        % Do nothing.
    case 'gravcapillary' % Capillary limit with gravity
        % Do nothing.
    otherwise
        error(['Method ''' method ''' not recognized.']);
end

end



