function solver = setUpGeothermalSolverAD(model)

    solver = AMGCLSolverAD('block_size', 2);
    solver.maxIterations = 100;
    solver.tolerance = 1e-4;
    solver.setPreconditioner('relaxation');
    solver.setRelaxation('ilu0');
    solver.setSolver('gmres');
    
%     solver.amgcl_setup.block_size = 2;
    nc = model.G.cells.num;
    solver.keepNumber = 2*nc;
    solver.reduceToCell = false;
    order = [1:nc; (1:nc) + nc];
    order = order(:);
%     
%     o = getCellMajorReordering(nc, 2);
%     
    [solver.equationOrdering, solver.variableOrdering] = deal(order);
    
    
end