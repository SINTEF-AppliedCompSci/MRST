function G=porosityadjustment(G)

nummatcell=G.Matrix.cells.num;
NNClist=G.nnc.cells;
fporevollist=G.nnc.porevol;
mporo=G.Matrix.rock.poro;
mvol=G.Matrix.cells.volumes;
count=0;
modcells=[];
for i=1:nummatcell
    connectedfracind=find(NNClist(:,1)==i);
    if isempty(connectedfracind), continue; end
    count=count+1;
    modcells(end+1)=i;
    matporevol=mporo(i)*mvol(i);
    fracporevol=sum(fporevollist(connectedfracind));
    totPV=fracporevol+matporevol;
    G.rock.poro(i)=totPV/mvol(i);
end

disp([num2str(count),' matrix cell porosities were modified.']);
% disp(modcells);




end