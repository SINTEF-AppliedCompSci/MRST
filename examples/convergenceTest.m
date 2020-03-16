mrstModule add mpfa mpsaw vemmech

Nd       = 2;
nref     = 2;
kappa    = 1;
alpha    = 1;
gridtype = 1;
eta      = 1e-8;

doVem = false;

output = mpsaPaperConvergenceFunc(Nd, nref, kappa, alpha, gridtype, eta, ...
                                  'doVem', doVem);

dex = output.dex; 
dnum = output.dnum;
deL2 = output.deL2;

fprintf('relative L2 error exact vs mpsa: %g\n', deL2(end));
if doVem
    fprintf('relative L2 error exact vs vem: %g\n', deVEM(end));
end


%% Print convergence rates

fprintf('Convergence rate for MPSA\n');
log2(deL2(1 : end - 1)./deL2(2 : end))
if doVem
    fprintf('Convergence rate for VEM\n');
    log2(deVEM(1 : end - 1)./deVEM(2 : end))
end

return
    
%% 
if opt.doVem
    n = 3; 
else 
    n = 2;
end

figure
set(gcf, 'numbertitle', 'off', 'name', 'DISPLACEMENT')

subplot(n, 2, 1)
title('x-mpsa')
plotCellData(G, dnum(: , 1), 'edgecolor', 'none'), colorbar

subplot(n, 2, 2)
title('y-mpsa')
plotCellData(G, dnum(: , 2), 'edgecolor', 'none'), colorbar

subplot(n, 2, 3)
title('x-exact')
plotCellData(G, dex(: , 1), 'edgecolor', 'none'), colorbar

subplot(n, 2, 4)
title('y-exact')
plotCellData(G, dex(: , 2), 'edgecolor', 'none'), colorbar

if opt.doVem
    subplot(n, 2, 5)
    title('x-vem')
    plotNodeData(G, uVEM(: , 1), 'edgecolor', 'none'), colorbar

    subplot(n, 2, 6)
    title('y-vem')
    plotNodeData(G, uVEM(: , 2), 'edgecolor', 'none'), colorbar
end
