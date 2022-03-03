function model = Fractional(model)
    muW = model.experiment.fluid.muW;
    muO = model.experiment.fluid.muNW;
    s = model.satfun.s;
    kr = model.satfun.kr;
    for i = 1 : length(s)
        krs = kr(s(i));
        krw = krs(1); kro = krs(2);
        M  = (krw / kro) * (muO / muW);
        fw(i) = 1 / (1 + 1 / M);
    end
    model.satfun.fw = fw';
end