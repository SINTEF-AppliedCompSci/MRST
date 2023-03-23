function index=findneighbours(G,cellind)
% Find the 26 neighbours around a cell
assert(isfield(G,'cartDims'));
N=getNeighbourship(G);


% first generation
index=any(N==cellind,2);
N1=N(index,:);
N1cellind=N1(~ismember(N1,cellind));

% second generation
N2cellind=zeros(36,1);
for i=1:length(N1cellind)
    cellindi=N1cellind(i);
    index=any(N==cellindi,2);
    Ni=N(index,:);
    Nicellind=Ni(~ismember(Ni,[cellindi,cellind]));
    N2cellind(6*(i-1)+(1:length(Nicellind)))=Nicellind;  
end
N2cellind=N2cellind(N2cellind>0);
bin=accumarray(N2cellind(:),1);
N2cellind=find(bin>1);
% disp(N2cellind);

% third generation
N3cellind=zeros(36,1);
for i=1:length(N2cellind)
    cellindi=N2cellind(i);
    index=any(N==cellindi,2);
    Ni=N(index,:);
    Nicellind=Ni(~ismember(Ni,[cellindi,cellind,N1cellind']));
    N3cellind(6*(i-1)+(1:length(Nicellind)))=Nicellind;  
end
N3cellind=N3cellind(N3cellind>0);
bin=accumarray(N3cellind(:),1);
N3cellind=find(bin>1);
% disp(N3cellind);


index=sort([N1cellind;N2cellind;N3cellind]);
end