for i = 1:1

    [G, rock, fluid, deck, state] = setupEGG('realization', i);
    
    figure; plotToolbar(G, rock);
end
 
 %%
 grd = readGridBoxArray(grd, fid, kw, nc);