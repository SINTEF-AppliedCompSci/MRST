G = cartGrid([30, 40]);

di = [5; 10; 10; 5];
dj = [5; 10; 10; 10; 5];

p = partitionTensor(G, di, dj);
clf;
plotCellData(G, mod(p, 13));

%%
p = partitionUniformPadded(G, [3, 4]);
clf;
plotCellData(G, mod(p, 13));
