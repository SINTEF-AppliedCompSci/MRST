function G = ascData2Grid(tdata,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
opt=struct('cartDims',tdata.grid_size-1,'nz',1,'FH',[],'Hdata',[]);
opt=merge_options(opt,varargin{:});
origo=tdata.grid_space([1,3])';
Lx=tdata.grid_space(2)-tdata.grid_space(1);
Ly=tdata.grid_space(4)-tdata.grid_space(3);
G=cartGrid([opt.cartDims,opt.nz],[Lx,Ly,1]);
G.nodes.coords(:,1:2)=bsxfun(@plus,G.nodes.coords(:,1:2),origo);

if(isempty(opt.FH))
    assert(~isempty(opt.Hdata))
    Hdata=opt.Hdata;
    % Hdata can ave nan
    x=linspace(Hdata.grid_space(1),Hdata.grid_space(2),Hdata.grid_size(1));
    %xx=(x(1:end-1)+x(2:end))/2;
    y=linspace(Hdata.grid_space(3),Hdata.grid_space(4),Hdata.grid_size(2));
    %yy=(y(1:end-1)+y(2:end))/2;
    [X,Y]=meshgrid(x,y);
    ind=isfinite(Hdata.val{3}');
    Xtmp=Hdata.val{1}';
    assert(all(X(ind)==Xtmp(ind)));
    Ytmp=Hdata.val{2}';
    assert(all(Y(ind)==Ytmp(ind)));
    H=Hdata.val{3}';
    fixed_value=max(H(ind));
    H(~ind)=max(fixed_value);
    FH=@(x,y) interp2(X,Y,H,x,y);
else
    fixed_value=100;
    FH=opt.FH;
end

% define top and bottum inerpolant
Ft=@(x,y) interp2(tdata.val{1}',tdata.val{2}',tdata.val{3}',x,y);
%FH=@(x,y) interp2(Hdata.val{1}',Hdata.val{2}',Hdata.val{3}',x,y);


zt=Ft(G.nodes.coords(:,1),G.nodes.coords(:,2));
assert(all(isfinite(zt)));
Hb=FH(G.nodes.coords(:,1),G.nodes.coords(:,2));
%Hb=Hb_int.Values;
%assert(all(isfinite(Hb)));
Hb(~isfinite(Hb))=fixed_value;

G.nodes.coords(:,3)=zt+G.nodes.coords(:,3).*Hb;

end

