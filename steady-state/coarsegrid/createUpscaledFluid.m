function fluid = createUpscaledFluid(f, updata, partition)
% Creates an upscaled fluid from given data.

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

% Create region
ntsat = max(partition);
if ntsat == 1
    reg.SATNUM = [];
    reg.SATINX = ':';
else
    reg.SATNUM = partition;
    reg.SATINX = cellfun(@(x)find(x==reg.SATNUM), num2cell(1:ntsat), ...
        'UniformOutput', false);
end


%--------------------------------------------------------------------------
% Two phase data
%--------------------------------------------------------------------------

if numel(updata(1).krW) > 1 || numel(updata(1).krO) > 1
    warning('SteadyState:MultidimRelperm', ...
        'Multidim relperm. Using only first dimension!');
end
s2c = @(fn) arrayfun(@(x) x.(fn){1}, updata, 'UniformOutput', false);


%--------------------------------------------------------------------------
% Relative permeability

TW = s2c('krW');
fluid.krW = @(sw, varargin) prop(TW,   sw, reg, varargin{:});
TO = s2c('krO');
fluid.krO = @(so, varargin) prop(TO, 1-so, reg, varargin{:});
fluid.krOW = fluid.krO;

%--------------------------------------------------------------------------
% Invert frac flow function

pref = 200*barsa; % viscosity reference pressure % TODO how to choose?
muW = f.muW(pref, 'cellInx', 1);
if isfield(f,'muO')
    muO = f.muO(pref, 'cellInx', 1);
else
    muO = f.BOxmuO(pref, 'cellInx', 1) / f.BO(pref, 'cellInx', 1);
end
% Fractional flow = (krW/muW) / ( (krW/muW) + (krO/muO) )
T = cell(size(TW));
for i=1:numel(TW)
    ff   = (TW{i}(:,2)/muW) ./ (TW{i}(:,2)/muW + TO{i}(:,2)/muO);
    T{i} = [ff, TW{i}(:,1)];
    % In case the fractional flow is not monotone, we skip the points which
    % do not increase in value.
    T{i} = T{i}([true; diff(T{i}(:,1))>0], :);
end
fluid.fracFlowInv = @(ff, varargin) prop(T, ff, reg, varargin{:});

%--------------------------------------------------------------------------
% Capillary pressure

if isfield(updata, 'pcOW')
    T = arrayfun(@(x) x.pcOW, updata, 'UniformOutput', false);
    fluid.pcOW = @(sw, varargin) prop(T, sw, reg, varargin{:});
    
    % Invert pcOW curve
    for i=1:numel(T)
        T{i} = [T{i}(:,2) T{i}(:,1)];
        if T{i}(1,1)>T{i}(2,1)
            T{i} = flipud(T{i}); % flip if reversed order
        end
    end
    fluid.pcOWInv = @(pc, varargin) prop(T, pc, reg, varargin{:});
end


% Copy remaining properties over to upscaled fluid structure
fns = fieldnames(f);
for i = 1:numel(fns)
    if ~isfield(fluid, fns{i}) % check existence
        fluid.(fns{i}) = f.(fns{i});
    end
end


end


%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------

% General frunction
function v = prop(T, sw, reg, varargin)
satInx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
T = extendTab(T);
v = interpReg(T, sw, satInx);  
end

% 
% % Water Relperm
% function v = krW(T, sw, reg, varargin)
% satInx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
% T = extendTab(T);
% v = interpReg(T, sw, satInx);  
% end
% 
% % Oil Relperm
% function v = krO(T, sw, reg, varargin)
% satInx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
% T = extendTab(T);
% v = interpReg(T, sw, satInx); 
% end
% 
% % Capillary Pressure
% function v = pcOW(T, sw, reg, varargin)
% satInx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
% T = extendTab(T);
% v = interpReg(T, sw, satInx);
% end



