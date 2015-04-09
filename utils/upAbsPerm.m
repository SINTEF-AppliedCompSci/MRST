function [KU, states] = upAbsPerm(G, varargin)
opt = struct(...
    'dp',         1*barsa, ...
    'Gp',         [], ...
    'bcp',        [], ...
    'areas',      [], ...
    'lengths',    [], ...
    'forcesym',   false, ...
    'absval',     false, ...
    'states',     [], ...
    'dims',       1:3, ...
    'bcFaces',    [], ...
    'psolver',    [], ...
    'returnFlux', false, ...
    'V',          [] ...
    );
opt = merge_options(opt, varargin{:});

blockInfo.lengths


%opt = checkOptions(G, opt);

% Run simulation(s) to find fluxes
V = findFluxes(G, rock, opt);

% Perform upscaling
KU = performUpscaling(V, opt);

end

% -------------------------------------------------------------------------


function V = findFluxes(B, opt)

dims  = opt.dims;
ndims = length(dims);

doReturnStates = nargout > 1;

% Initial state
state0 = initResSol(B.G, 100*barsa, 1);

% Flux matrix / vector
V = nan(ndims,1);

% Choose solver if not given
fluidPure  = initSingleFluid('mu' ,1, 'rho', 1);
if B.periodic
    opt.psolver = @(state, rock, bcp) incompTPFA(state, B.G, ....
        computeTransGp(B.G.parent, B.G, rock), fluidPure, 'bcp', bcp);
else
    opt.psolver = @(state0, rock, bc) incompTPFA(state0, B.G, ...
        computeTrans(G, rock), fluidPure, 'bc', bc);
end

% Apply pressure drop in each dimension
for i = 1:ndims
    
    % Set pressure drop
    if B.periodic
        opt.bcp.value(:) = 0;
        opt.bcp.value(opt.bcp.tags == dims(i)) = opt.dp(i);
        bc = opt.bcp;
    else
        bc = addBC([], opt.bcFaces{dims(i),1}, 'pressure', opt.dp(i));
        bc = addBC(bc, opt.bcFaces{dims(i),2}, 'pressure', 0);
    end
    
    % Solve
    warning('off','mrst:periodic_bc');
    state1 = opt.psolver(state0, rock, bc);
    warning('on','mrst:periodic_bc');
    
    if doReturnStates
        states{i} = state1;
    end
    
    if B.periodic
        % Store flux in j-direction caused by pressure drop in d-direction
        for j = 1:ndims
            faces = opt.bcp.face(opt.bcp.tags==dims(j));
            sign  = opt.bcp.sign(opt.bcp.tags==dims(j));
            V(j,i) = sum(state1.flux(faces, 1).*sign) / opt.areas(j);
        end
    else
        faces = opt.bcFaces{dims(i), 2};
        sign  = ones(numel(faces), 1);
        sign(G.faces.neighbors(faces,1)==0) = -1;
        V(i) = sum(state1.flux(faces, 1).*sign) / opt.areas(i);
    end
    
end


end


function KU = performUpscaling(V, opt)

isPeriodic = ~isempty(opt.Gp);
dims  = opt.dims;
ndims = length(dims);
L  = opt.lengths;
dp = opt.dp;


% Compute upscaled permeability vector
if isPeriodic
    Pinv = diag(L(:)./dp(:));
    KU   = - V*Pinv;
else
    KU   = V.*(L(:)./dp(:));
end

% Check for negative elements
if opt.absval
    if any(any(KU<0))
        KU = abs(KU);
    end
end

% Enforce symmetry of KU if requested
if opt.forcesym && isPeriodic
    if ndims > 1
        a = (KU(1,2) + KU(2,1))/2; KU(1,2) = a; KU(2,1) = a;
    end
    if ndims > 2
        a = (KU(1,3) + KU(3,1))/2; KU(1,3) = a; KU(3,1) = a;
        a = (KU(2,3) + KU(3,2))/2; KU(2,3) = a; KU(3,2) = a;
    end
end

end




function opt = checkOptions(G, opt)

ndims = length(opt.dims);

isPeriodic = ~isempty(opt.Gp);
if isPeriodic
    assert(~isempty(opt.bcp), 'bcp must be given for periodic grids.');
end

% Compute grid lengths if not given
if isempty(opt.lengths)
    L = max(G.faces.centroids) - min(G.faces.centroids);
    L = L(opt.dims);
else
    L = opt.lengths;
    assert(length(L) == ndims, ['Option ''lengths'' must same length '...
        'as ''dims''.']);
end
opt.lengths = L;

opt.dp = opt.dp .* ones(ndims,1);

assert(max(opt.dims) <= G.griddim, 'Dimensions does not match grid.');
assert(length(opt.dp) == 1 || length(opt.dp) == ndims, ...
    'Parameter dp must be scalar or vector of same length as dimensions.');

% Gravity not supported
assert(norm(gravity) == 0, 'Upscaling with gravity not supported.');

% Find sides of fixed grid if not given
if ~isPeriodic && isempty(opt.bcFaces)
    opt.bcFaces = getBoundaryFaces(G);
end

% Choose solver if not given
if isempty(opt.psolver)
    fluidPure  = initSingleFluid('mu' ,1, 'rho', 1);
    if isPeriodic
        opt.psolver = @(state, rock, bcp) incompTPFA(state, opt.Gp, ....
            computeTransGp(G, opt.Gp, rock), fluidPure, 'bcp', bcp);
    else
        opt.psolver = @(state0, rock, bc) incompTPFA(state0, G, ...
            computeTrans(G, rock), fluidPure, 'bc', bc);
    end
end

% Compute side areas if not given
if isempty(opt.areas)
    opt.areas = zeros(ndims, 1);
    for i = 1:ndims
        if isPeriodic
            opt.areas(i) = sum( opt.Gp.faces.areas( ...
                opt.bcp.face(opt.bcp.tags==opt.dims(i))) );
        else
            opt.areas(i) = sum( G.faces.areas(opt.bcFaces{opt.dims(i),2}) );
        end
    end
else
    assert(length(opt.areas) == ndims, ...
        'Option ''areas'' must have same length as dims.');
end

% Allocate space for fluxes
if isempty(opt.V)
    if isPeriodic
        opt.V = nan(ndims,ndims); % matrix
    else
        opt.V = nan(ndims,1); % vector
    end
end

end



