% Testing the flooding
clear all
close all

nx = 10;
physdims = [1, 1, 1];
G0 = cartGrid([nx, nx, nx], physdims);
G = computeGeometry(G0);
n = [1, 1, 0];

%% Cut along diag
[G2, ix2, g2] = sliceGrid2(G, physdims/2, 'normal', n);
start = 0.1*physdims;
m1 = flood(G2, start, ix2.new.faces);
start = 0.9*physdims;
m2 = flood(G2, start, ix2.new.faces);
assert(all(m1+m2==1))
figure
plotCellData(G2, m1+2*m2);
view(3)

%% Cut almost along diag
offset = 1e-2;
[G2, ix2, g2] = sliceGrid2(G, physdims/2+offset, 'normal', n);
start = 0.1*physdims;
m1 = flood(G2, start, ix2.new.faces);
start = 0.9*physdims;
m2 = flood(G2, start, ix2.new.faces);
assert(all(m1+m2==1))
figure
plotCellData(G2, m1+2*m2);
view(3)

%% Two parallel cuts, three domains
cuts = [physdims/2; physdims/2+offset];
[G2, ix2, g2] = sliceGrid2(G, cuts, 'normal', n);
start = 0.1*physdims;
m1 = flood(G2, start, ix2.new.faces);
start = 0.9*physdims;
m2 = flood(G2, start, ix2.new.faces);
start = mean(cuts);
m3 = flood(G2, start, ix2.new.faces);
assert(all(m1+m2+m3==1));
figure
plotCellData(G2, m1+2*m2+3*m3);
view(3)
