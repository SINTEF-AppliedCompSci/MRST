% Test marking cut domains from sliceGrid

clear all
close all

legacy = verLessThan('matlab', '8.6');

nx = 5;
physdims = [1, 1, 1];
G0 = cartGrid([nx, nx, nx], physdims);
G0 = computeGeometry(G0);
n = [1, 1, 0];

%% Single cut
[G, gix] = sliceGrid(G0, physdims/2, 'normal', n);
start = 0.1*physdims;
if legacy
    m = markCutGrids(G, gix.new.faces, 'legacy', true, 'start', [0.1,0.1,0.1]);
else
    m = markCutGrids(G, gix.new.faces);
end
figure
plotCellData(G, m);
view(3)

%% Two parallel cuts, three domains
offset = 0.1;
cuts = [physdims/2; physdims/2+offset];
[G, gix] = sliceGrid(G0, cuts, 'normal', n);
if legacy
    x0 = [0.1, 0.1, 0.1];
    m1 = markCutGrids(G, gix.new.faces, 'legacy', true, 'start', 1-x0);
    m2 = markCutGrids(G, gix.new.faces, 'legacy', true, 'start', x0);
    m = m1+2*m2;
else
    m = markCutGrids(G, gix.new.faces);
end
figure
plotCellData(G, m);
view(3)
