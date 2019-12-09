function model = addMultiscaleSolver(model, CG, varargin)
    s = getSmootherFunction('type', 'ilu', 'iterations', 1);
    
    mssolver = MultiscaleVolumeSolverAD(CG, 'getSmoother', s, varargin{:});
    msmodel = MultiscalePressureModel(model.pressureModel.G,...
                                      model.pressureModel.rock,...
                                      model.pressureModel.fluid,...
                                      model.pressureModel, mssolver);
    model.pressureModel = msmodel;
end