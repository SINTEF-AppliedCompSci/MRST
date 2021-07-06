function setup = updateSetupFromScaledParameters(setup, params, p)
nparam = cellfun(@(x)x.nParam, params);
p      = mat2cell(p, nparam, 1);
for k = 1:numel(params)
    pu  = params{k}.unscale(p{k});
    setup = params{k}.setParameterValue(setup, pu);
end
end