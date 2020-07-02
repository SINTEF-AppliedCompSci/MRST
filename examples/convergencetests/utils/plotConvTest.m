function plotConvTest(output, params)
    
    deL2 = output.deL2;
    
    Nd = params.Nd;
    nref = params.nref;
    kappa = params.kappa;
    alpha = params.alpha;
    gridtype = params.gridtype;
    eta = params.eta;

    log2N = (1 : nref)';
    plot(log2N, log2(deL2), '*-', 'linewidth', 4);
    ax = gca;
    ax.YAxis.FontSize = 18;
    ax.XAxis.FontSize = 18;
    ylabel('log2(error)', 'fontsize', 18);
    xlabel('-log2(1/N)', 'fontsize', 18);
    caseTitle = setCaseTitle(Nd, gridtype, eta, kappa, alpha);
    title(caseTitle);    
end


function casetitle = setCaseTitle(Nd, gridtype, eta, kappa, alpha)
    dimstr = sprintf('%dD', Nd);
    switch gridtype
      case 1
        gridtypestr = 'Cartesian grid';
      case 2
        gridtypestr = 'Triangular grid, 90 degree angles';
      case 3
        gridtypestr = 'Equilateral triangles (alt)';
      case 4
        gridtypestr = 'Equilateral triangles';
      case 5
        gridtypestr = 'Tetrahedrals';
    end
    
    casetitle = sprintf('%s - %s, \\kappa = %0.3g, \\alpha = %0.3g, \\eta = %0.3g', dimstr, ...
                        gridtypestr, kappa, alpha, eta);
        
end


