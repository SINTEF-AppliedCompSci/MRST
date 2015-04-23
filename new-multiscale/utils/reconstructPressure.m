function [p, solvetime] = reconstructPressure(CG, pressure, A, rhs)
    D = formReconstructionMatrix(A, CG.partition);
    
    tmp1 = warning('query','MATLAB:nearlySingularMatrix');
    tmp2 = warning('query','MATLAB:singularMatrix');
    warning('off','MATLAB:nearlySingularMatrix')
    warning('off','MATLAB:singularMatrix')
    t = tic();
    p = mldivide(D, rhs - (A - D)*pressure);
    solvetime = toc(t);
    warning(tmp1.state,'MATLAB:nearlySingularMatrix')
    warning(tmp2.state,'MATLAB:singularMatrix')
end
