function s = setupSimCompVe(G, rock, varargin);
  rock_tmp=rock;
    rock_tmp.perm=rock.perm.*G.cells.H;
    T = computeTrans(G,rock_tmp);
    % take out internal trans
    N  = double(G.faces.neighbors);
    intInx = (prod(N,2)~=0);
    cf = G.cells.faces(:,1); 
    nf = G.faces.num;
    T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
    T_all=T;
    T  = T(intInx);
    %    
    pv = poreVolume(G,rock);
    pv = pv.*G.cells.H;    
    s = setupSimComp(G, rock, 'porv', pv,'trans',T);
    s.T_all=T_all;
end