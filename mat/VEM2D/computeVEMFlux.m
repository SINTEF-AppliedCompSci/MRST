function state = computeVEMFlux(G, rock, state, bc)

T = computeTrans(G, rock);

ii = 1:G.faces.num;



ii = G.cells.faces(:,1);
jj = 1:numel(G.cells.faces(:,1));

P = sparse(ii, jj, 1);
T = 1./(P*(1./T));

bf = boundaryFaces(G);
f = true(G.faces.num,1);
f(bf) =false;


pDiff = zeros(G.faces.num,1);
pDiff(f) = state.cellPressure(G.faces.neighbors(f,1)) ...
                        - state.cellPressure(G.faces.neighbors(f,2));
for i = 1:numel(bf)
    
    f = bc.face(i);
    v = bc.value(i);
    
    if strcmp(bc.type(i), 'pressure')
        if G.faces.neighbors(f,1) == 0
            pDiff(f) = v - state.cellPressure(G.faces.neighbors(f,2));
        else
            pDiff(f) = state.cellPressure(G.faces.neighbors(f,1))-v;
        end
    else
        T(f) = 1;
        pDiff(f) = v*(-1).^(G.faces.neighbors(f,1) == 0);
    end 
end

flux = T.*pDiff;

state.flux = flux;

end