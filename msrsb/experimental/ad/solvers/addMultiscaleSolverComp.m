function model = addMultiscaleSolverComp(model, CG, varargin)
    opt = struct('MultiscaleModel', 'compositional');
    [opt, extra] = merge_options(opt, varargin{:});
    s = getSmootherFunction('type', 'ilu', 'iterations', 1);
    mssolver = MultiscaleVolumeSolverAD(CG, 'getSmoother', s, 'maxIterations', 5,...
                                            'tolerance', 1e-3, extra{:});
    
    if isa(model.pressureModel, 'WrapperModel')
        parent = model.pressureModel.parentModel;
    else
        parent = model.pressureModel;
    end
        
    msmodel = CompositionalMSPressureModel(parent.G,...
                                      parent.rock,...
                                      parent.fluid,...
                                      parent, mssolver);
    model.pressureModel = msmodel;
    model.pressureModel.multiscaleSolver.resetBasis = false;
    model.pressureModel.multiscaleSolver.updateBasis = true;
    model.pressureModel.multiscaleSolver.updateInterval = 10;
end