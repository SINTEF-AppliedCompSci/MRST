function model = addMultiscaleSolverComp(model, CG, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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
