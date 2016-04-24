function [Gl,T] = makeNNCextruded(G,Gl,F,fracture,flayers)
%
nfl = numel(flayers);

% frac-matrix connections
G = frac_matrix_nnc(G,F,fracture);
% Remove unnecessary connections
if any(G.nnc.CI==0)
    G.nnc.cells = G.nnc.cells(G.nnc.CI~=0,:);
    G.nnc.type(G.nnc.CI==0,:) = [];
    G.nnc.CI = G.nnc.CI(G.nnc.CI~=0);
end
nm = G.cells.num;
findex2D = []; findex3D = [];
for i = 1:numel(F)
    findex2D = [findex2D,F(i).cells.start]; %#ok
    findex3D = [findex3D, ...
        Gl.FracGrid.(['Frac',num2str(i)]).cells.start]; %#ok
end
findex2D(end+1) = findex2D(end) + F(i).cells.num;
findex3D(end+1) = findex3D(end) + Gl.FracGrid.(['Frac',num2str(i)]).cells.num;
Gl.nnc.cells = []; Gl.nnc.CI = repmat(G.nnc.CI,nfl,1);
Gl.nnc.type = repmat(G.nnc.type,nfl,1);
for i = 1:size(G.nnc.cells,1)
    mat = G.nnc.cells(i,1);
    frac = G.nnc.cells(i,2);
    plane = find(frac>=findex2D(1:end-1) & frac<findex2D(2:end));
    ind2D = frac - findex2D(plane);
    ind3D = findex3D(plane) + ind2D;
    temp = [mat + (flayers'-1)*nm, ind3D + ((1:nfl)'-1)*F(plane).cells.num];
    Gl.nnc.cells = [Gl.nnc.cells;temp];
end
    
Gl = assembleGlobalGrid(Gl); 
Gl.fractureAperture = fracture.aperture;

if isfield(Gl.rock,'poro'), pv = poreVolume(Gl,Gl.rock);
else pv = Gl.cells.volumes; end
w1 = pv(Gl.nnc.cells(:,1))./Gl.rock.perm(Gl.nnc.cells(:,1));
w2 = pv(Gl.nnc.cells(:,2))./Gl.rock.perm(Gl.nnc.cells(:,2));
wt = pv(Gl.nnc.cells(:,1))+pv(Gl.nnc.cells(:,2));
% w1 = 1./Gl.rock.perm(Gl.nnc.cells(:,1));
% w2 = 1./Gl.rock.perm(Gl.nnc.cells(:,2));
% wt = 1;
Gl.nnc.T = Gl.nnc.CI.*(wt./(w1+w2)); clear wt w1 w2

Gl = frac_frac_nnc_extruded(G, Gl, F, fracture, flayers);
Gl = rmfield(Gl,{'numLayers','layerSize','type'});
% compute transmissibilities
T = computeTrans(Gl, Gl.rock);
%-------------------------------------------------------------------------%
% computeTrans returns 2 transmissibilities for each internal face and one
% transmissibility fo each external face. Size of transmissibility array is
% same as G.cells.faces if opt.usetrans = false 
%-------------------------------------------------------------------------%
cf = Gl.cells.faces(:,1);
nf = Gl.faces.num;
T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
T = [T;Gl.nnc.T];

Gl.faces.tag = zeros(Gl.faces.num,1);
x = ismember(Gl.faces.neighbors,Gl.nnc.cells,'rows');
Gl.faces.tag(x) = 1;
return