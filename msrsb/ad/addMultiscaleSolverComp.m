function model = addMultiscaleSolverComp(model, CG, varargin)
    opt = struct('MultiscaleModel', 'compositional');
    [opt, extra] = merge_options(opt, varargin{:});
    s = getSmootherFunction('type', 'ilu', 'iterations', 1);
    mssolver = MultiscaleVolumeSolverAD(CG, 'getSmoother', s, 'maxIterations', 5,...
                                            'tolerance', 1e-3, extra{:});
    pmodel = model.pressureModel;
    if isa(model.pressureModel, 'WrapperModel')
        parent = pmodel.parentModel;
    else
        parent = pmodel;
    end
    msmodel = CompositionalMSPressureModel(parent.G,...
                                      parent.rock,...
                                      parent.fluid,...
                                      pmodel, mssolver);
    model.pressureModel = msmodel;
    model.pressureModel.multiscaleSolver.resetBasis = false;
    model.pressureModel.multiscaleSolver.updateBasis = true;
    model.pressureModel.multiscaleSolver.updateInterval = 10;
end