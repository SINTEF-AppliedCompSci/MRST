function fluid = initADIFluidOWPolymer(varargin)
% Make a structure representing a three-component fluid (water, oil,
% polymer) and their properties (relative permeabilities, densities,
% viscosities, etc.). This is might be a helpful function for examples and
% testing.
%
% SYNOPSIS:
%   fluid = initADIFluidOWPolymer()
%   fluid = initADIFluidOWPolymer('pn1', 'pv1', ...)
%
% PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.
%             The supported options are as follows:
%   regnum  - Fluid index.
%             The grid is divided into N regions, where each region has
%             different fluid properties. This parameter is a vector of
%             length G.cells.num, such that cell i belongs to fluid area
%             reginx(i). If empty, entire grid has the same properties.
%
% FLUID PROPERTY PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.
%   Each property value can be given in of one of the following four forms:
%     1) Constant in entire grid.
%        Input form: P=k, where k is a scalar.
%        The property value will be set to k in all cells.
%     2) Constant in each region.
%        Input form: P=v, where v is a vector of size N-by-1, where N is
%        the number of regions,
%        The property value will be set to v(i) in all cells of region i.
%     3) Piecewise linear function, same function in entire domain.
%        Input form: P=T, where T is matrix of size M-by-2. The linear
%        function is defined by the M points (xi,yi), where xi=T(i,1) and
%        yi=T(i,2)).
%     4) Piecewise linear function, different function in each region.
%        Input form: P={T1,...,TN}, where Ti is matrix of size M-by-2,
%        defining the linear function for region i. N is the number of
%        regions.
%     5) Piecewise linear function of two variables (x,v), same function in
%        entire domain.
%        Input form: P=T, where T is structure with the following fields:
%           T.key  - The v values.
%           T.pos  - Maps vi values to the data field, such that
%                       T.data(T.pos(i):(T.pos(i+1)-1), :)
%                    gives the x and y values corresponding to vi=T.key(i).
%           T.data - A matrix of size M-by-2, where xj=T.data(j,1) and
%                    yj=T.data(j,2).
%
%   Note that if a property value is empty, the property is not added to
%   the fluid structure.
%
%   The supported properties are as follows:
%
%   Water/Oil properties:
%   muW     - Water viscosity (Pa*s), function of water pressure
%   muO     - Oil viscosity (Pa*s), function of oil pressure
%   rhoW    - Water density (kg/m^3)
%   rhoO    - Oil density (kg/m^3)
%   BW      - Water formation volume factor (FVF), function of water
%             pressure
%   BO      - Oil formation volume factor (FVF), function of oil pressure
%   krW     - Water relative permeability, function of water saturation
%   krO     - Oil relative permeability, function of oil saturation
%   pvMultR - Pressure-dependent pore volume multiplier, function of oil
%             pressure.
%
%   Polymer properties:
%   cmax    - Maximum polymer concentration (kg/m^3)
%   ads     - Adsorption isoterm (kg/kg), function of polymer
%             concentration (kg/m^3)
%   adsMax  - Maximum allowed adsorbed polymer (kg/kg)
%   muWMult - Water viscosity multiplier (kg/m^3), function of polymer
%             concentration (kg/m^3). muWMult is defined such that
%               muM(c) = muW.*muWMult(c)
%             where muM(c) is the viscosity of a fully mixed polymer
%             solution.
%   dps     - Dead pore space
%   rrf     - Residual resistance factor
%   rhoR    - Mass density of the rock formation (kg/m^3)
%   mixPar  - Todd-Longstaff mixing parameter omega
%
% RETURNS:
%   fluid   - ADI fluid structure with the given properties.
%

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

require ad-core ad-props

opt = struct(...
    'satnum',   [],  ...
    'muW',      1, ...
    'muO',      1, ...
    'rhoW',     1, ...
    'rhoO',     1, ...
    'BW',       1, ...
    'BO',       1, ...
    'isIncomp', false, ...
    'krW',      1, ...
    'krO',      1, ...
    'pcOW',     [], ...
    'pvMultR',  [], ...
    'ads',      [], ...
    'adsMax',   [], ...
    'adsInx',   [], ...
    'muWMult',  [], ...
    'dps',      [], ...
    'rrf',      [], ...
    'rhoR',     [], ...
    'mixPar',   [], ...
    'cmax',     [], ...
    'swir',     [], ...
    'sor',      [], ...
    'fracFlowInv', true, ...
    'pcOWInv',     true, ...
    'polymer',     true ... % set to false to ignore polymer fields
    );

% If an odd number of arguments are given, we assume the first argument is
% a structure with options. We convert the structure to varargin form and
% pass it to mergeOptions.
if mod(nargin,2)==1
    p = varargin{1};
    assert(isstruct(p), 'A single argument must be a structure.');
    vararginp = struct2args(p);
    opt = mergeOptions(opt, vararginp{:});
    varargin = varargin(2:end); % The remaining arguments
end

opt = mergeOptions(opt, varargin{:});



% Handle regions
regnum = opt.satnum;
reginx = getReginx(regnum); % cell array with indecies for each domain
if isempty(regnum)
    nreg = 1;
else
    nreg   = max(regnum); % number of different fluid areas
end

% Set empty fluid
fluid  = [];



% Oil / Water properties --------------------------------------------------

phases = {'W', 'O'};
for j = 1:numel(phases)
    n     = phases{j};
    addProperty(['mu'  n],     opt.(['mu'  n]), true);
    addProperty(['mu'  n 'S'], opt.(['mu'  n]), true);
    addProperty(['rho' n],     opt.(['rho' n]), false);
    addProperty(['rho' n 'S'], opt.(['rho' n]), false);
    addProperty(['B'   n],     opt.(['B'   n]), true);
    fluid.(['b' n]) = @(varargin) 1./fluid.(['B' n])(varargin{:});
    addProperty(['kr'  n],     opt.(['kr'  n]), true);
    
    % TODO: Not sure how the following are used
    %addProperty(['krO' n],     opt.(['kr'  n]), true);
end
fluid.rhoGS = 0;

% Add relperm field
fluid.relPerm = @(sW, varargin) relPerm(sW, varargin{:});


addProperty('pcOW', opt.pcOW, true);
%fluid.relPerm = @(sw, varargin) relPerm(krW, krO, sw, varargin{:});
addProperty('pvMultR', opt.pvMultR,     true);

% Saturation limits
%addProperty('swir', opt.swir, false);
%addProperty('sor', opt.sor, false);
if ~isempty(opt.swir)
    assert(numel(opt.swir)==1 || numel(opt.swir)==nreg);
    fluid.swir = opt.swir;
end
if ~isempty(opt.sor)
    assert(numel(opt.sor)==1 || numel(opt.sor)==nreg);
    fluid.sor = opt.sor;
end

% Satnum
fluid.satnum = opt.satnum;

% Incompressible
if opt.isIncomp
    fluid.isIncomp = true;
end

% Inverse fractional flow curve
if opt.fracFlowInv
    %    % Fractional flow = (krW/muW) / ( (krW/muW) + (krO/muO) )
    %    krW = opt.krW; krO = opt.krO;
    %    muW = opt.muW; muO = opt.muO;
    %    assert(isscalar(muW) && isscalar(muO), 'Assumption failed.');
    %    wascell = iscell(krW);
    %    if ~wascell
    %       krW = {krW}; krO = {krO};
    %    end
    %    T = cell(1, numel(krW));
    %    for j = 1:numel(krW)
    %       T{j} = [( krW{j}(:,2)/muW )./( (krW{j}(:,2)/muW) + ...
    %          (flipud(krO{j}(:,2))/muO) ), krW{j}(:,1)];
    %    end
    %    T = extendTab(T);
    %    if ~wascell
    %       T = T{1};
    %    end
    %    addProperty('fracFlowInv', T, true);
    
    fluid.fracFlowInv = @(ff, varargin) fracFlowInv(ff, opt, ...
        regnum, reginx, varargin{:});
    
    if opt.polymer
        fluid.fracFlowInvPoly = @(ff, c, varargin) fracFlowInvPoly(ff, ...
            c, opt, regnum, reginx, varargin{:});
    end
    
end

% Inverse capillary pressure curve
if ~isempty(opt.pcOW) && opt.pcOWInv
    pcOW = opt.pcOW;
    wascell = iscell(pcOW);
    if ~wascell
        pcOW = {pcOW};
    end
    T = cell(1, numel(pcOW));
    for j = 1:numel(pcOW)
        T{j} = [pcOW{j}(:,2), pcOW{j}(:,1)];
        if T{j}(1,1) > T{j}(2,1)
            T{j} = flipud(T{j});
        end
    end
    T = extendTab(T);
    if ~wascell
        T = T{1};
    end
    addProperty('pcOWInv', T, true);
end



% Polymer properties ------------------------------------------------------

if opt.polymer
    addProperty('ads',     opt.('ads'),     true);
    addProperty('adsMax',  opt.('adsMax'),  false);
    addProperty('adsInx',  opt.('adsInx'),  false);
    addProperty('muWMult', opt.('muWMult'), true);
    addProperty('dps',     opt.('dps'),     false);
    addProperty('rrf',     opt.('rrf'),     false);
    addProperty('rhoR',    opt.('rhoR'),    false);
    addProperty('mixPar',  opt.('mixPar'),  false);
    addProperty('cmax',    opt.('cmax'),    false);
end



% NESTED FUNCTIONS --------------------------------------------------------

function addProperty(name, P, isFunc)
    % Add a property to the fluid structure
    %   name   - Name of property in fluid structure
    %   P      - The property input
    %   isFunc - Boolean value indicating if the property is a function taking
    %            an input parameter.

    if isempty(P)
        % Do not add the property

    elseif iscell(P) && ischar(P{1})

        % Special case
        key = P{1};
        val = P{2};
        addSpecialProperty(name, key, val, isFunc);

    elseif isnumeric(P) && isscalar(P)

        % Case 1) Constant in entire grid.

        % Check input
        assert(~isnan(P), 'Constant is NaN');

        if isFunc
            %fluid.(name) = @(x, varargin) P;
            fluid.(name) = @(x, varargin) constantProperty(...
                repmat(P,nreg,1), x, regnum, {':'}, varargin{:});
        else
            fluid.(name) = P;
        end

    elseif isnumeric(P) && size(P,2)==1

        % Case 2) Constant in each region.

        % Check input
        assert(~any(isnan(P)), 'Some values are NaN');
        assert(size(P,1)==nreg, ...
            ['Property ''' name ''' is of incorrect length.']);

        if isFunc
            fluid.(name) = @(x, varargin) ...
                constantProperty(P, [], regnum, reginx, varargin{:});
        else
            %          fluid.(name) = @(varargin) ...
            %             constantProperty(P, [], regnum, reginx, varargin{:});
            fluid.(name) = constantProperty(P, [], regnum, reginx);
        end

    elseif isnumeric(P) && size(P,2)==2

        % Case 3) Piecewise linear function, same function in entire domain.

        % Check input
        assert(~any(any(isnan(P))), 'Some values are NaN');
        assert(isFunc, ...
            ['Property ''' name ''' does not take an input parameter and ' ...
            'can therefore not be given as a piecewise linear function.']);
        assert(numel(unique(P(:,1))) == size(P,1), ...
            ['Property ''' name ''': x-values must be unique.']);


        P = extendTab(P);
        fluid.(name) = @(x, varargin) piecewiseLinearCurve(...
            P, x, [], [], varargin{:});

    elseif iscell(P) && numel(P)==nreg && isnumeric(P{1})

        % Case 4) Piecewise linear function, different func in each region.

        % Check input
        assert(~any(cellfun(@(x) any(any(isnan(x))), P)), ...
            'Some values are NaN');
        assert(isFunc, ...
            ['Property ''' name ''' does not take an input parameter and ' ...
            'can therefore not be given as piecewise linear functions.']);
        for i = 1:nreg
            assert(numel(unique(P{i}(:,1))) == size(P{i},1), ...
                ['Property ''' name ''': x-values must be unique (region ' ...
                num2str(i) ').']);
            P{i} = extendTab(P{i});
        end

        fluid.(name) = @(x, varargin) piecewiseLinearCurve(P, x, ...
            regnum, reginx,  varargin{:});

    elseif isstruct(P)

        % Case 5) Piecewise linear function of two variables (x,v), same
        % function in entire domain.
        assert(isFunc, ...
            ['Property ''' name ''' does not take any input parameters ' ...
            'and can therefore not be given as a function.']);
        fluid.(name) = @(x, v, varargin) piecewiseLinear2DCurve(...
            P, x, v, [], [], varargin{:});

    elseif iscell(P) && numel(P)==nreg && struct(P{1})

        % Case 6) Piecewise linear function of two variables (x,v). Different
        % function in each domain.

        assert(isFunc, ...
            ['Property ''' name ''' does not take any input parameters ' ...
            'and can therefore not be given as a function.']);
        fluid.(name) = @(x, v, varargin) piecewiseLinear2DCurve(...
            P, x, v, regnum, reginx, varargin{:});

    else

        % Form of the input is not supported. Throw error.
        error(['Input to property ''' name ''' not given correctly.']);

    end

end



function addSpecialProperty(name, key, props, isFunc)
    % Add a property to the fluid structure
    %   name   - Name of property in fluid structure
    %   P      - The property input
    %   isFunc - Boolean value indicating if the property is a function taking
    %            an input parameter.

    switch key
        case 'pvcdo'
            % This is the Formation Volume Factor for oil as described in the
            % ECLIPSE manual in the PVCDO keyword.
            assert(strcmpi(name, 'BO'), ['The PVCDO special property '...
                'may only be used for ''BO''.']);
            assert(numel(props)>3-1, 'PVCDO parameters not given correctly.');
            assert(numel(props)<3+1, 'PVCDO only supporting three parameters.');
            assert(isFunc, ...
                ['Property ''' name ''' does not take any input parameters ' ...
                'and can therefore not be given as a function.']);
            fprintf('Adding special property PVCDO.\n');

            % Get values and setup function handle
            pref   = convertFrom(props(1), barsa); % Reference pressure
            fvfref = convertFrom(props(2), 1); % BO at ref pressure
            comp   = convertFrom(props(3), 1/barsa); % Compressibility
            fluid.BO = @(p,varargin) fvfref.*exp(-comp.*(p-pref));

        case 'pvtw'
            % This is the Formation Volume Factor for water as described in the
            % ECLIPSE manual in the PVTW keyword.
            assert(strcmpi(name, 'BW'), ['The PVTW special property '...
                'may only be used for ''BW''.']);
            assert(numel(props)>3-1, 'PVTW parameters not given correctly.');
            assert(numel(props)<3+1, 'PVTW only supporting three parameters.');
            assert(isFunc, ...
                ['Property ''' name ''' does not take any input parameters ' ...
                'and can therefore not be given as a function.']);
            fprintf('Adding special property PVTW.\n');

            % Get values and setup function handle
            pref   = convertFrom(props(1), barsa); % Reference pressure
            fvfref = convertFrom(props(2), 1); % BW at ref pressure
            comp   = convertFrom(props(3), 1/barsa); % Compressibility
            X = @(p) comp.*(p - pref);
            fluid.BW = @(p,varargin) fvfref./( 1 + X(p) + X(p).^2/2 );

        case 'rock'
            % This is the Rock Compressibility formulation as described in the
            % ECLIPSE manual in the ROCK keyword. This is a multiplier for the
            % pore volume at the reference pressure.
            %
            % We assume there is only one PVT region, i.e. that NTPVT=1.
            %
            % Default values: 1.0132 barsa, 4.934E-5 1/barsa
            assert(strcmpi(name, 'pvMultR'), ['The ROCK special property '...
                'may only be used for ''pvMultR''.']);
            assert(numel(props)>2-1, 'ROCK input not given correctly.');
            assert(numel(props)<2+1, 'Only supporting two parameters in ROCK.');
            assert(isFunc, ...
                ['Property ''' name ''' does not take any input parameters ' ...
                'and can therefore not be given as a function.']);
            fprintf('Adding special property ROCK.\n');

            % Get values and setup function handle
            pref   = convertFrom(props(1), barsa); % Reference pressure
            crref  = convertFrom(props(2), 1/barsa); % Rock compr at ref pres
            X = @(p) crref.*(p - pref);
            fluid.pvMultR = @(p,varargin) 1 + X(p) + X(p).^2/2;

        otherwise
            error(['Special property key ''' key ''' unknown.']);
    end

end



function [krW, krO] = relPerm(sW, varargin)
    krW = fluid.krW(sW, varargin{:});
    krO = fluid.krO(1 - sW, varargin{:});
end

end


%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------


function prm = mergeOptions(prm, varargin)
% This function is an edited version of the file merge_options from the
% MATLAB Reservoir Simulation Toolbox (MRST).
%
% The difference is that this function accepts a new value of a different
% class than the default value.

if nargin > 1,
    if mod(numel(varargin), 2) == 0 && ...
            all(iscellstr(varargin(1 : 2 : end))),
        st = dbstack(1);
        try
            caller = st(1).name;
        catch  %#ok
            caller = 'BASE';
        end
        ofn = fieldnames(prm);
        nfn = varargin(1 : 2 : end);
        nfv = varargin(2 : 2 : end);
        
        for n = 1 : numel(nfn),
            ix = find(strcmpi(nfn{n}, ofn));
            
            if ~isempty(ix),
                prm.(ofn{ix}) = nfv{n};
            else
                warning([caller, ':Option:Unsupported'], ...
                    ['Option `', nfn{n}, ''' is not supported']);
            end
        end
    else
        error(msgid('Input:Huh'), ...
            'Huh? Did you remember to unpack VARARGIN?!?');
    end
end
end



function reginx = getReginx(regnum)
% Given regnum, get reginx.

if isempty(regnum)
    reginx = ':';
else
    nregs = max(regnum);
    reginx = cellfun(@(x)find(x==regnum), num2cell(1:nregs), ...
        'UniformOutput', false);
end

end





function yi = constantProperty(Y, xi, regnum, reginx, varargin)
% Y vector of length nreg
% N number of yi values to extract

nreg = numel(reginx);

if nreg > 0 && ischar(reginx{1}) && strcmp(reginx{1}, ':'),
    % Special case denoting entire domain in single region.
    
    if ~isempty(xi) % return yi same length as xi
        if isa(xi, 'ADI')
            N = numel(xi.val);
        else
            N = numel(xi);
        end
    else
        N  = numel(regnum); % all cells in grid
    end
    yi = repmat(Y(1), N, 1);
    
elseif nreg > 0
    
    rinx = getRegMap(xi, regnum, reginx, varargin{:});
    
    % Count the length of yi
    N = 0;
    for k = 1:nreg
        N = N +  numel(rinx{k});
    end
    yi = zeros(N,1);
    
    for k = 1:nreg,
        yi(rinx{k}) = repmat(Y(k), numel(rinx{k}), 1);
    end
    
end


end



function yi = piecewiseLinearCurve(T, xi, regnum, reginx, varargin)
% Create a piecewise linear curve.
%
% SYNOPSIS:
%   val   = piecewiseLinearCurve(T, xi, regnum, reginx)
%   val   = piecewiseLinearCurve(T, xi, regnum, reginx, 'pn1', 'pv1', ...)
%
% PARAMETERS:
%   T       - cell array of points for the curve
%   xi      -
%   regnum  -
%   reginx  -
%
% OPTIONAL PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.
%             The supported options are as follows:
%   cellInx -
%
% RETURNS:
%   yi      - Interpolated values for each cell in the grid, or for each
%             cell in cellInx, if given.
%

if isa(xi, 'ADI')
    N = numel(xi.val);
else
    N = size(xi,1);
end

if isempty(reginx)
    % Use T in entire domain
    yi  = interpReg({T}, xi, {(1:N)'} );
else
    %inx = getRegMapL(N, regnum, reginx, varargin{:});
    
    % hack for scalar input for entire domain
    if isnumeric(xi) && isscalar(xi) && isempty(varargin)
        xi = repmat(xi, numel(regnum), 1);
    end
    
    inx = getRegMap(xi, regnum, reginx, varargin{:});
    
    % if isscalar(xi) && N > 1
    %    % Allow scalar input if cellInx is given
    %    xi = repmat(xi, N, 1);
    % end
    yi  = interpReg(T, xi, inx);
end

end


function yi = piecewiseLinear2DCurve(T, xi, vi, regnum, reginx, varargin)
% Create a 2D piecewise linear curve.
%
% SYNOPSIS:
%   yi = piecewiseLinearCurve(T, xi, yi, regnum, reginx)
%   ti = piecewiseLinearCurve(T, xi, yi, regnum, reginx, 'pn1', 'pv1', ...)
%
% PARAMETERS:
%   T       - cell array of points for the curve
%   xi      -
%   vi      -
%   regnum  - Cell array
%   reginx  - Vector
%
% OPTIONAL PARAMETERS:
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.
%             The supported options are as follows:
%   cellInx -
%
% RETURNS:
%   yi      - Interpolated values for each cell in the grid, or for each
%             cell in cellInx, if given.
%

if isa(xi, 'ADI')
    N = numel(xi.val);
else
    N = size(xi,1);
end

flag = false(N,1);

if isempty(reginx)
    % Use T in entire domain
    yi = interpRegPVT({T}, xi, vi, flag, {(1:N)'});
else
    %inx = getRegMapL(N, regnum, reginx, varargin{:});
    
    % hack for scalar input for entire domain
    if isnumeric(xi) && isscalar(xi) && isempty(varargin)
        xi = repmat(xi, numel(reginx), 1);
    end
    if isnumeric(vi) && isscalar(vi) && isempty(varargin)
        vi = repmat(vi, numel(reginx), 1);
    end
    
    inx = getRegMap(xi, regnum, reginx, varargin{:});
    
    yi = interpRegPVT(T, xi, vi, flag, inx);
end


end




% Inverse fractional flow curve
function ffInv = fracFlowInv(ff, opt, regnum, reginx, varargin)

% Fractional flow = (krW/muW) / ( (krW/muW) + (krO/muO) )
krW = opt.krW; krO = opt.krO;
muW = opt.muW; muO = opt.muO;
assert(isscalar(muW) && isscalar(muO), 'Assumption failed.');

if iscell(krW)
    
    % Case 4) Piecewise linear function, different func in each region.
    T = cell(1, numel(krW));
    for j = 1:numel(krW)
        T{j} = [( krW{j}(:,2)/muW )./( (krW{j}(:,2)/muW) + ...
            (flipud(krO{j}(:,2))/muO) ), krW{j}(:,1)];
    end
    T = extendTab(T);
    ffInv = piecewiseLinearCurve(T, ff, regnum, reginx, varargin{:});
    
else
    
    % Case 3) Piecewise linear function, same function in entire domain.
    T = [( krW(:,2)/muW )./( (krW(:,2)/muW) + ...
        (flipud(krO(:,2))/muO) ), krW(:,1)];
    T = extendTab(T);
    ffInv = piecewiseLinearCurve(T, ff, [], [], varargin{:});
    
end

end

% Fractional flow, where krWtilde = krW / Rk
function ffInv = fracFlowInvPoly(ff, c, opt, regnum, reginx, ...
    varargin)

% Fractional flow = (krW/muW) / ( (krW/muW) + (krO/muO) )
assert(isscalar(c));

% Get asorption
ads = opt.ads;
ads = extendTab(ads);
if iscell(ads)
    adsVal = nan(numel(ads),1);
    for i = 1:numel(ads)
        adsVal(i) = interp1(ads{i}(:,1), ads{i}(:,2), c);
    end
else
    adsVal = interp1(ads(:,1), ads(:,2), c);
end

% Get max adsorption
adsMax = opt.adsMax;

% Compute Rk (scalar or vector)
Rk = 1 + (opt.rrf(:) - 1) .* (adsVal./adsMax(:));

% Scale the relative permeability (in each region) with 1/Rk.
krW = opt.krW;
if iscell(krW)
    assert(numel(Rk)==numel(krW));
    for i = 1:numel(krW)
        krW{i}(:,2) = krW{i}(:,2) ./ Rk(i);
    end
else
    assert(isscalar(Rk));
    krW(:,2) = krW(:,2) ./ Rk;
end

opt.krW = krW;

ffInv = fracFlowInv(ff, opt, regnum, reginx, varargin{:});

end
