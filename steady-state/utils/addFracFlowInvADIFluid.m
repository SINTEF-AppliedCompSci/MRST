function fluid = addFracFlowInvADIFluid(fluid, deck, varargin)
% Add a fluid function for computing the inverse of the water fractional
% flow curve. This may be called later by using
%   sW = fluid.fracFlowInv(ff)
% The input deck must have the properties SWOF, PVTW and PVCDO.
% 
% Fractional flow = (krW/muW) / ( (krW/muW) + (krO/muO) )
% 


T = createFracFlowTablesFromDeck(deck);
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

