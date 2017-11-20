function plot2DPoroExample(model, states, opt)
% Helper function for runAll2Dcases
    figure()
    plotToolbar(model.G, states);
    title(sprintf('fluid model: %s, method: %s', opt.fluid_model, opt.method));
end