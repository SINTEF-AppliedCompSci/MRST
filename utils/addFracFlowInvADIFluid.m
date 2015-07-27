function fluid = addFracFlowInvADIFluid(fluid, deck, varargin)
% Add a fluid function for computing the inverse of the water fractional
% flow curve. This may be called later by using
%   sW = fluid.fracFlowInv(ff)
% The input deck must have the properties SWOF, PVTW and PVCDO.
% 
% Fractional flow = (krW/muW) / ( (krW/muW) + (krO/muO) )
% 


assert(isfield(deck, 'PROPS'), 'Invalid deck. Missing ''PROPS'' field.');
assert(isfield(deck.PROPS, 'SWOF'), 'Missing property ''SWOF''.');
assert(isfield(deck.PROPS, 'PVTW'), 'Missing property ''PVTW''.');
%assert(isfield(deck.PROPS, 'PVCDO'), 'Missing property ''PVCDO''.');

% PVTW
% Columns: [REF.PRES., REF.FVF, COMPRESSIBILITY, REF.VISC., VISCOSIBILITY]
muW   = deck.PROPS.PVTW(4); % water viscosity
if deck.PROPS.PVTW(5) > 0
   warning(['There is nonzero water viscosibility. Only the reference '...
      'viscosity is used.']);
end

% PVCDO
if isfield(deck.PROPS, 'PVCDO')
   % Columns: [REF.PRES., REF.FVF, COMPRESSIBILITY, REF.VISC., VISCOSIBILITY]
   muO   = deck.PROPS.PVCDO(4); % oil viscosity
   if deck.PROPS.PVCDO(5) > 0
      warning(['There is nonzero water viscosibility. '...
          'Only the reference viscosity is used.']);
   end
elseif isfield(deck.PROPS, 'PVDO')
    muO = deck.PROPS.PVDO{1}(1,3); % oil viscosity
elseif isfield(fluid, 'muO')
   %muO = fluid.muO(100*barsa
   error('Have to be implemented')
end

% SWOF
% Columns: [Sw, Krw, Krow, Pcow]
% Fractional flow = (krW/muW) / ( (krW/muW) + (krO/muO) )
T  = cellfun(@(x) [( x(:,2)/muW )./( (x(:,2)/muW) + (x(:,3)/muO) ), ...
   x(:,1)], deck.PROPS.SWOF, 'UniformOutput', false);
T = extendTab(T);

% Add to fluid
reg   = handleRegions(deck, varargin{:});
fluid = assignFracFlowInv(fluid, T, reg);

end

% Remaining code is a slight modification of assignSWOF

function f = assignFracFlowInv(f, T, reg)
f.fracFlowInv = @(ff, varargin)FracFlowInv(ff, T, reg, varargin{:});
end


function v = FracFlowInv(ff, T, reg, varargin)
satinx = getRegMap(ff, reg.SATNUM, reg.SATINX, varargin{:});
v = interpReg(T, ff, satinx);
end

