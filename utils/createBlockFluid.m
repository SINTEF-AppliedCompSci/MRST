function BFluid = createBlockFluid(fluid, cells, varargin)
% Extracting the fluid for the current coarse cell only.
opt = struct(...
    'cellInx',     false ...
    );
opt = merge_options(opt, varargin{:});

BFluid = fluid;

% Properties that take 'cellInx' name first
p = {'krW', 'krO', 'krOW', 'relPerm', 'pcOW', 'BW', 'bW', 'muW', ...
    'muO', 'muWMult', 'ads', 'pcOWInv', 'fracFlowInv'};
for i=1:numel(p)
    if isfield(fluid, p{i})
        BFluid.(p{i}) = @(x, varargin) fluid.(p{i})(x, ...
          'cellInx', subcells(cells, varargin{:}));
    end
end

% Properties that DO NOT normally take 'cellInx' name first
p = {'BO', 'bO', 'BOxmuO'};
for i=1:numel(p)
    if isfield(fluid, p{i})
        if opt.cellInx
            BFluid.(p{i}) = @(x, varargin) fluid.(p{i})(x, ...
                'cellInx', subcells(cells, varargin{:}));
        else
            BFluid.(p{i}) = @(x, varargin) fluid.(p{i})(x, ...
                subcells(cells, varargin{:}));
        end
    end
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


