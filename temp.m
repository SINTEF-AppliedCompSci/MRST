clear seltbl
seltbl.nodes = 5;
seltbl.num = 1;

seltbl = crossTable(seltbl, nodefacecoltbl, {'nodes'});

newnfctbl = addLocInd(nodefacecoltbl, 'nfcind');

fds = {'nodes', 'faces', 'coldim'};
seltbl = crossTable(seltbl, newnfctbl, fds);
ind = seltbl.nfcind;