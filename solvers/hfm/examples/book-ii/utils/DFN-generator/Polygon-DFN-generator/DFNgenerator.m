function [fracplanes,P32gen]=DFNgenerator(fracplanes,fracinput,physdim,exclzonemult,tol,varargin)
% NOTE: size inputs are radii. but length exponent is diameter
opt=struct('displayexpdomain',false,'displayDFN',false,'periodic',true,'circle',false);
opt=merge_options(opt,varargin{:});

opt.circle=fracinput.circle;

if opt.periodic
    [fracplanes,P32gen]=periodicDFNfracplanesgenerator(fracplanes,fracinput,physdim,exclzonemult,tol,...
        'displayDFN',opt.displayDFN,'circle',opt.circle);
else
    error('Work in progress...');
end

% [fracplanes.circle]=deal(opt.circle);


end