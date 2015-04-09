function fluid = createUpscaledFluid(f, updata, partition)
% Creates an upscaled fluid from given data.

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


% Two phase data ----------------------------------------------------------

s2c = @(fn) arrayfun(@(x) x.(fn), updata, 'UniformOutput', false);

T = s2c('krW');
fluid.krW = @(sw, varargin) krW(T,   sw, reg, varargin{:});
T = s2c('krO');
fluid.krO = @(so, varargin) krO(T, 1-so, reg, varargin{:});
fluid.krOW = fluid.krO;

if isfield(updata, 'pcOW')
    T = s2c('pcOW');
    fluid.pcOW = @(sw, varargin) pcOW(T, sw, reg, varargin{:});
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

function c = struct2CellArray(s, fn)
assert(isstruct(s) && isfield(s, fn));
c = cell(numel(s),1);

end


% Water Relperm
function v = krW(T, sw, reg, varargin)
satInx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
T = extendTab(T);
v = interpReg(T, sw, satInx);  
end

% Oil Relperm
function v = krO(T, sw, reg, varargin)
satInx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
T = extendTab(T);
v = interpReg(T, sw, satInx); 
end

% Capillary Pressure
function v = pcOW(T, sw, reg, varargin)
satInx = getRegMap(sw, reg.SATNUM, reg.SATINX, varargin{:});
T = extendTab(T);
v = interpReg(T, sw, satInx);
end



