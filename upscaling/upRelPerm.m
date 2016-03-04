function [updata, report] = upRelPerm(block, updata, ...
    method, varargin)
% 
opt = struct(...
    'nvalues',     20, ...
    'viscousmob',  false, ... % viscous upscaling: use total mob method
    'dims',        1:3, ...  % Dimensions to upscale
    'dp',          1*barsa, ...  % Pressure drop
    'savesat',     false, ... % save saturation distributions
    'absmethod',   'pressure' ... % one-phase upscaling method
    );
opt = merge_options(opt, varargin{:});

wantReport  = nargout > 1;
wantSatDist = opt.savesat && wantReport;
timeStart   = tic;

dims  = opt.dims;
ndims = length(dims);
nvals = opt.nvalues;
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

% Get input values depending on the method
values = getValues(block, updata, method, nvals);

% Loop over the input values
for iv = 1:nvals
    
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
        sW = directionDistribution(block, method, sW0, d, opt.dp);
        
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
        
        if strcmpi(method, 'viscous') && opt.viscousmob
            % For the viscous limit upscaling, we may use the total
            % mobility and only call the one phase upscaling once.
            
            % NOTE: We have experienced some issues with this method that
            % we have not solved. Upscaling of the SPE10 case got some
            % wrong results with this method.
            
            if wantReport
                report.viscousmob = true;
            end
            
            % TODO how to chose pressure?
            pref = 200*barsa; % viscosity reference pressure
            muW = fluid.muW(pref, 'cellInx', 1);
            if isfield(fluid,'muO')
                muO = fluid.muO(pref, 'cellInx', 1);
            else
                try
                    muO = fluid.BOxmuO(pref, 'cellInx', 1) / ...
                        fluid.BO(pref, 'cellInx', 1);
                catch
                    % Fallback to support old syntax
                    muO = fluid.BOxmuO(pref) / fluid.BO(pref);
                end
            end
            
            mobTot = kr{1}./muO + kr{2}./muW;
            rock_KmobT = rock;
            
            % If the permeability is anisotropic, we multiply each
            % direction with the total mobility
            rock_KmobT.perm = bsxfun(@times, rock.perm, mobTot);
            
            if all(rock_KmobT.perm == 0)
                % We will get no fluid motion
                KMobTU = 0;
            else
                if any(rock_KmobT.perm < 0)
                    error('Some pseudo perm values are negative!');
                end
                if any(any(rock_KmobT.perm == 0))
                    % To avoid a singular matrix when performing one
                    % phase upscaling, the zero permeabilities are set
                    % to something larger than zero.
                    ep = eps(mean(mean(...
                        rock_KmobT.perm(rock_KmobT.perm>0)) ));
                    rock_KmobT.perm(rock_KmobT.perm == 0) = ep*1e1;
                end
                
                % Perform one phase upscaling with the altered
                % permeability field
                block.rock = rock_KmobT;
                data = upAbsPerm(block, [], 'dims', d, 'dp', opt.dp, ...
                    'method', opt.absmethod);
                KMobTU = data.perm;
            end
            
            % For viscous limit, the value is fractional flow
            ffW = val;
            krwat = muW * ffW * KMobTU / Kup(id); % water
            kroil = muO * (1 - ffW) * KMobTU / Kup(id); % oil
            krO{id}(iv,:) = [sWup, kroil];
            krW{id}(iv,:) = [sWup, krwat];
            
        else
            
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
        
        end
        
    end % End of dimension loop
    
end % End of input value loop

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
    report.nvalues = nvals;
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
        assert(all(isfield(block.fluid, {'swir','sor'})), ...
            ['Fluid structure needs fields ''swir'', ''sor'' and ' ...
            '''satnum''.']);
        
        swir   = block.fluid.swir;
        sor    = block.fluid.sor;
        satnum = block.fluid.satnum;
        pvTot  = sum(block.pv);
        sWUmin = sum(swir(satnum).*block.pv)/pvTot;
        sWUmax = sum((1-sor(satnum)).*block.pv)/pvTot;
        
        % Equal spacing of the upscaled saturations
        values = linspace(sWUmin, sWUmax, nvals)';
        
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
        
    case 'viscous' % Vuscous limit
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


function sW = directionDistribution(block, method, sW, d, dp)
% Update saturation for the current direction

% Get saturation values
switch method
    case 'flow' % Flow based simulation
        
        %error('Need update') % TODO: We need to get d here.
        % TODO: Also, should we update the saturation in each direction?
        
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
            % Cells neighbouring the outer faces
            nc = @(i) sum(block.G.faces.neighbors(block.faces{d}{i},:),2);
            
            sWf = sW(nc(1));
            bc = addBC([], block.faces{d}{1}, ...
                'pressure', dp, 'sat', [sWf 1-sWf]);
            
            sWf = sW(nc(2));
            bc = addBC(bc, block.faces{d}{2}, ...
                'pressure',  0, 'sat', [sWf 1-sWf]);
            
            Gp  = [];
            bcp = [];
        end
        
        % Initialize solution
        state0      = initResSol(G, 100*barsa, [sW, 1-sW]);
        state0.flux = zeros(G.faces.num, 2);
        
        % Solve steady state flow
        state1 = simulateToSteadyStateADI_new(G, block.rock, ...
            block.fluid, sW, 'bc', bc, 'Gp', Gp, 'bcp', bcp);
        
        % Extract saturation
        sW1 = state1.s(:,1);
        
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





