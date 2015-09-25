% assess sensitivity of grid


res1 = testSleipnerSensFUN('refineLevel',-4);
res2 = testSleipnerSensFUN('refineLevel',-4,'addPerturbation','true','pertAmp',2);


figure
subplot(1,2,1)
title('unperturbed top')
plotCellData(res1.Gt, res1.dobj_dz, 'EdgeColor','none')
axis equal tight
colorbar
hold on
line(res1.plumes{end}.outline(:,1), res1.plumes{end}.outline(:,2), 'LineWidth',2,'Color','r')
subplot(1,2,2)
title('perturbed top')
plotCellData(res2.Gt, res2.dobj_dz, 'EdgeColor','none')
axis equal tight
colorbar
hold on
line(res2.plumes{end}.outline(:,1), res2.plumes{end}.outline(:,2), 'LineWidth',2,'Color','r')

hfig = gcf;