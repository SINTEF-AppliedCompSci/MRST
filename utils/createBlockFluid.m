function BFluid = createBlockFluid(fluid, cells)
% Extracting the fluid for the current coarse cell only.

BFluid = fluid;

% Properties that take 'cellInx' name first
p = {'krW', 'krO', 'krOW', 'relPerm', 'pcOW', 'BW', 'bW', 'muW', ...
    'muO', 'muWMult', 'ads', 'pcOWInv', 'fracFlowInv'};
for i=1:numel(p)
    BFluid.(p{i}) = @(x, varargin) fluid.(p{i})(x, ...
      'cellInx', subcells(cells, varargin{:}));
end

% Properties that DO NOT take 'cellInx' name first
p = {'BO', 'bO', 'BOxmuO'};
for i=1:numel(p)
    BFluid.(p{i}) = @(x, varargin) fluid.(p{i})(x, ...
        subcells(cells, varargin{:}));
end

end


%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------

function c = subcells(cells, varargin)
% This function is used to extract the correct cells from a fluid property.

opt = struct('cellInx', []);
opt = merge_options(opt, varargin{:});
if isempty(opt.cellInx)
   c = cells;
else
   c = cells(opt.cellInx);
end

end


