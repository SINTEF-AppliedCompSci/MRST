function plotNornePoroExample(model, states, opt)
% Helper function for runAll2Dcases
    figure()
    plotToolbar(model.G, states, 'outline', true);
    view([7, 40]);
    title(sprintf('fluid model: %s, method: %s', opt.fluid_model, opt.method));
end