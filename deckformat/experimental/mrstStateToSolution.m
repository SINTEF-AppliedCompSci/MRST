function SOLUTION = mrstStateToSolution(G,state0)

ind=G.cells.indexMap;

press = zeros(prod(G.cartDims),1);
press(ind) = state0.pressure/barsa;

s=zeros(prod(G.cartDims),3);
s(ind,:)=state0.s;

rs=zeros(prod(G.cartDims),1);
rs(ind)=state0.rs;

rv=zeros(prod(G.cartDims),1);
rv(ind)=state0.rv;
SOLUTION=struct('PRESSURE',press,'SWAT',s(:,1),'SGAS',s(:,3),'RS',rs,'RV',rv)
end