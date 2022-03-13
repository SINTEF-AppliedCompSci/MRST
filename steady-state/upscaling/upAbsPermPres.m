function [updata, report] = upAbsPermPres(block, updata, varargin)
%Undocumented Utility Function

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
    'dp',         1*barsa, ...
    'dims',       1:3, ...
    'psolver',    'tpfa', ...
    'fulltensor', false ... % Return full tensor or just the diagonal
    );
opt = merge_options(opt, varargin{:});

if nargin==1
    updata = [];
end

wantReport = nargout > 1;
timeStart = tic;

% Handle input
dims  = opt.dims; % Dimensions to upscale
ndims = length(dims);
dp    = opt.dp .* ones(ndims,1);
G     = block.G;
rock  = block.rock;
isPeriodic = block.periodic;

% Initial state
state0 = initResSol(G, 100*barsa, 1);

% Allocate space for fluxes
if isPeriodic
    V = nan(ndims, ndims); % matrix
else
    V = nan(ndims,1); % vector
end

% Setup solver
fluidPure  = initSingleFluid('mu' ,1, 'rho', 1);

switch opt.psolver
    % Mimetic is always consistent, but TPFA seems to be faster.
	case 'mimetic'
        if isPeriodic
            bcp = block.bcp;
            S = computeMimeticIPGp(G.parent, G, rock);
            psolver = @(state0, bcp) incompMimetic(state0, G, S, fluidPure, ...
                'bcp', bcp);
        else
            S = computeMimeticIP(G, rock);
            psolver = @(state0, bc) incompMimetic(state0, G, S, fluidPure, ...
                'bc', bc);
        end
    case 'tpfa'
        if isPeriodic
            bcp = block.bcp;
            T = computeTransGp(G.parent, G, rock);
            psolver = @(state0, bcp) incompTPFA(state0, G, ....
                T, fluidPure, 'bcp', bcp);
        else
            T = computeTrans(G, rock);
            psolver = @(state0, bc) incompTPFA(state0, block.G, ...
                T, fluidPure, 'bc', bc);
        end
    otherwise
        error('Pressure solver type ''%s'' unknown.', opt.psolver);
end

% Loop over dimensions, apply pressure drop and compute fluxes
for i = 1:ndims
    
    % Set boundary conditions
    if isPeriodic
        bcp.value(:) = 0;
        bcp.value(bcp.tags == dims(i)) = dp(i);
        bc = bcp;
    else
        bc = addBC([], block.faces{dims(i)}{1}, 'pressure', dp(i) );
        bc = addBC(bc, block.faces{dims(i)}{2}, 'pressure', 0 );
    end
    
    % Solve pressure equation
    warning('off','mrst:periodic_bc');
    warning('off','all');
    state1 = psolver(state0, bc);
    warning('on','all');
    warning('on','mrst:periodic_bc');
    
    if isPeriodic
        % Store flux in j-direction caused by pressure drop in d-direction
        for j = 1:ndims
            faces = bcp.face(bcp.tags==dims(j));
            sign  = bcp.sign(bcp.tags==dims(j));
            V(j,i) = sum(state1.flux(faces, 1).*sign) / ...
                block.areas(dims(j),2);
        end
    else
        faces = block.faces{dims(i)}{2};
        sign  = ones(numel(faces), 1);
        sign(G.faces.neighbors(faces,1)==0) = -1;
        V(i) = sum(state1.flux(faces, 1).*sign) / block.areas(dims(i),2);
    end
    
end

Kup = [];

% Compute upscaled permeability
L = block.lengths(dims);
if isPeriodic
    Pinv = diag(L(:)./dp(:)); % matrix
    Kup  = - V*Pinv;
    if ~opt.fulltensor
        Kup  = diag(Kup); % The diagonal is extracted
    end
else
    Kup   = V.*(L(:)./dp(:)); % vector
    if opt.fulltensor
        Kup = diag(Kup); % Create a tensor matrix
    end
end
if ~opt.fulltensor
    % If only a vector is returned, it is returned as a row-vector
    Kup = Kup(:)';
end

updata.perm = Kup;

totalTime = toc(timeStart);
if wantReport
    report.psolver  = opt.psolver;
    report.periodic = isPeriodic;
    report.dims     = dims;
    report.dp       = dp;
    report.time     = totalTime;
end

end
