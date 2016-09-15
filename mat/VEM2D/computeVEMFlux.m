function state = computeVEMFlux(G, rock, state, bc)

T = computeMultiPointTrans(G, rock);
T  = T.T;
bf  = any(G.faces.neighbors==0, 2);
I1 = [(1:G.cells.num)'; G.cells.num + find(bf)];
p = [state.cellPressure; state.facePressure(bf)];
state.flux = cellFlux2faceFlux(G, T(:, I1) * p);

end