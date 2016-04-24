function [Ab, bb] = extractSubBlock(CG, state, fluid, p_ms, nearWell, W)
G = CG.parent;
state.pressure = p_ms;
xr = state;
[mu, ~] = fluid.properties(xr);
s         = fluid.saturation(xr);
kr        = fluid.relperm(s, xr);

mob    = bsxfun(@rdivide, kr, mu);
totmob = sum(mob, 2);
T = computeTrans(G,G.rock);
cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
T = T.*totmob(cellNo);
p_ms_face =  ...
    accumarray(G.cells.faces(:,1), p_ms(cellNo).*T, [G.faces.num,1])./ ...
    accumarray(G.cells.faces(:,1), T, [G.faces.num,1]);
[subG,gc,gf,~] = extractSubgrid(G,nearWell);
for i = 1:numel([W.cells])
    wcells = find(gc == W(i).cells);
    W(i).cells = wcells;
end
boundary = find(any(subG.faces.neighbors==0,2));
gbound = any(G.faces.neighbors==0,2);
% in_boundary = boundary(~ismember(gf(boundary),find(gbound)));
% out_boundary = boundary(ismember(gf(boundary),find(gbound)));
% b_in_global = gf(in_boundary);
% b_out_global = gf(out_boundary);
% 
% bcfaces = [in_boundary;out_boundary]; %in_boundary;%
% bcval_in = zeros(size(in_boundary));
% bcval_out = zeros(size(out_boundary));
% 
% for i = 1:numel(in_boundary)
%     fnbrs = G.faces.neighbors(b_in_global(i),:);
%     bcval_in(i) = p_ms(fnbrs(~ismember(fnbrs,nearWell)));
% end
% for i = 1:numel(out_boundary)
%     fnbrs = G.faces.neighbors(b_out_global(i),:);
%     bcval_out(i) = p_ms(fnbrs(fnbrs~=0));
% end
% bcval = [bcval_in;bcval_out]; %bcval_in; %
% bc = addBC([],bcfaces,'pressure',bcval);
boundary = boundary(~ismember(gf(boundary),find(gbound)));
bc = addBC([],boundary,'pressure',p_ms_face(gf(boundary)));

subG.rock.perm = G.rock.perm(nearWell);
subG.rock.poro = G.rock.poro(nearWell);
sT = computeTrans(subG,subG.rock);
subState = initResSol (subG, state.pressure(gc), state.s(gc));
[Ab,bb] = getSystemIncompTPFA(subState, subG, sT, fluid, 'bc', bc, 'Wells', W);
end