function u = well2control(W, varargin)
% Produce control vector u from target-wells
% If scaling is supplied, u is scaled s.t. 0<=u<=1
% according to scaling.boxLims
opt = struct('targets', (1:numel(W))', ...
             'scaling', []);
opt  = merge_options(opt, varargin{:});
vals = vertcat(W(opt.targets).val);
if ~isempty(opt.scaling)
    lims = opt.scaling.boxLims;
    u = (vals - lims(:,1))./(lims(:,2)-lims(:,1));
    if any(or(u<-eps, u>1+eps))
        warning('Some well-values lie outside given box-constraints');
    end
else
    u = vals;
    warning('No scaling was given, setting control-values equal to target well-values');
end
end
