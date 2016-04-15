%%
files=dir('data/*.asc');
for k=1:numel(files)
    data=readSeismicAsc(fullfile('data',files(k).name));
    data.file=files(k);
    %{
    data.F=@(vec) interp2(data.val{1}',data.val{2}',data.val{3}',vec(:,1),vec(:,2));
    data.Fc=@(X,Y) interp2(data.val{1}',data.val{2}',data.val{3}',X,Y)
    %}
    data.F=@(vec) interp2(data.X,data.Y,data.Z,vec(:,1),vec(:,2));
    data.Fc=@(X,Y) interp2(data.X,data.Y,data.Z,X,Y)
    datas{k}=data;
end
%%
smodelsnames={'IEAGHGmodel','ORIGINALmodel'}
smodels={};
for i=1:numel(smodelsnames)
[ G, Gt, rock, rock2D ] = makeSleipnerModelGrid('modelName',smodelsnames{i});
Fttmp=scatteredInterpolant(Gt.cells.centroids(:,1),Gt.cells.centroids(:,2),Gt.cells.z);
Ft=@(vec) Fttmp(vec(:,1),vec(:,2));
Fc=@(X,Y) reshape(Fttmp(X(:),Y(:)),size(X));
 smodels{i}=struct('G',G,'Gt',Gt,'Ft',Ft,'Fc',Fc);
end
%% utsira model
[grdecls datasets petroinfo] = getAtlasGrid('utsirafm');
meta=datasets{2}.meta;
data=datasets{2}.data;
Gu = processAAIGrid(meta, data, false, [1 1]);
Ftu=scatteredInterpolant(Gu.nodes.coords(:,1),Gu.nodes.coords(:,2),Gu.nodes.coords(:,3));
%{
dims = [meta.ncols, meta.nrows];
h = meta.cellsize;
xx=h*[0:(dims(1)-1)];
yy=h*[0:(dims(2)-1)];
[xx,yy]=ndgrid(xx,yy);
xx=xx+meta.xllcorner;
yy=yy+meta.yllcorner;
zz=(abs(data'));
Ftu=@(X,Y) interp2(xx,yy,zz,X,Y);
%}
%%
%%
dpad=[1e3,1e3]*0;
res=[200,200];
%res=[800,1200];
n=1;
[X,Y]=meshgrid(linspace(datas{n}.grid_space(1)-dpad(1),datas{n}.grid_space(2)+dpad(1),res(1)),...
    linspace(datas{n}.grid_space(3)-dpad(2),datas{n}.grid_space(4)+dpad(2),res(2)));

%%
fprintf('\n')
for i=1:numel(datas)
   fprintf('%i \t %s\n',i,datas{i}.file.name); 
end
%%
% 8-6 is approximatly 7ms
% time conversion is approximatly constant (3*2)/(1*4) but varies in space
% 3-1 seems related to CO2 depth in ms for 2006
% 8-1 seems to be related but not the CO2 depth at 2008
% 7  is the co2deph in meters in 2008
% 5+6 is smoothly related to 7 and is some CO2 hight in ms
% what is relation ship between 8 and 6 ??

% not ok 5,6,7
%%
% 5 seams to be in ms with negative sign 
% what is relation between 1 and 6
%
%mesh(X,Y,Fttmp(X,Y))
mesh(X,Y,datas{2}.Fc(X,Y)-Ftc(X,Y))
%mesh(X,Y,datas{4}.Fc(X,Y)-Ftc(X,Y))
mesh(X,Y,-datas{8}.Fc(X,Y)+datas{1}.Fc(X,Y)+datas{5}.Fc(X,Y)*0+4)
%%
% 2008 Co2 depth
%mesh(X,Y,datas{5}.Fc(X,Y)+datas{6}.Fc(X,Y))% seems to be co2depth
clf,mesh(X,Y,(datas{5}.Fc(X,Y)+datas{6}.Fc(X,Y))./datas{7}.Fc(X,Y))% relation to depht??
xv=reshape(X,[],1);
yv=reshape(Y,[],1);
vv=reshape((datas{5}.Fc(X,Y)+datas{6}.Fc(X,Y))./datas{7}.Fc(X,Y),[],1);
ind=isfinite(vv);
ys=abs(min(yv(ind)));
xs=abs(min(xv(ind)));
dxs=max(xv(ind))-min(xv(ind));
dys=max(yv(ind))-min(yv(ind));
yv=(yv-ys)/dys;
xv=(xv-xs)/dxs;
a=[ones(sum(ind),1),yv(ind),yv(ind).^2,xv(ind)]\vv(ind);
hold on;
mesh(X,Y,a(1)+a(2)*(Y-ys)/dys+a(3)*((Y-ys).*(Y-ys))/dys.^2+a(4)*(X-xs)/dxs)
%mesh(X,Y,datas{5}.Fc(X,Y)+datas{6}.Fc(X,Y)-datas{7}.Fc(X,Y))
%plot(datas{5}.Fc(X,Y)+datas{6}.Fc(X,Y),datas{7}.Fc(X,Y),'*')
%plot(datas{5}.Fc(X,Y)+datas{6}.Fc(X,Y),datas{7}.Fc(X,Y),'*')
%max(max((datas{5}.Fc(X,Y)+datas{6}.Fc(X,Y))./datas{7}.Fc(X,Y)))
%%
plot(-datas{8}.Fc(X,Y)+datas{1}.Fc(X,Y)+datas{5}.Fc(X,Y),datas{5}.Fc(X,Y)+datas{6}.Fc(X,Y),'*')
% the depth conversion do not seem to change in time
mesh(X,Y,(datas{3}.Fc(X,Y).*datas{2}.Fc(X,Y))./(datas{4}.Fc(X,Y).*datas{1}.Fc(X,Y)))
mesh(X,Y,(datas{8}.Fc(X,Y).*datas{2}.Fc(X,Y))./(datas{4}.Fc(X,Y).*datas{6}.Fc(X,Y)))
mesh(X,Y,(datas{8}.Fc(X,Y)./datas{6}.Fc(X,Y)))
mesh(X,Y,(datas{8}.Fc(X,Y)./(datas{6}.Fc(X,Y)+0*datas{6}.Fc(X,Y))))
mesh(X,Y,(datas{8}.Fc(X,Y)-datas{6}.Fc(X,Y)));%&+0*datas{6}.Fc(X,Y))))
figure(1)
%mesh(X,Y,(datas{1}.Fc(X,Y)-datas{1}.Fc(X,Y))+(datas{6}.Fc(X,Y)+datas{5}.Fc(X,Y)+datas{1}.Fc(X,Y)-datas{8}.Fc(X,Y)));%&+0*datas{6}.Fc(X,Y))))
mesh(X,Y,-(datas{8}.Fc(X,Y)-datas{1}.Fc(X,Y))+7+0*datas{6}.Fc(X,Y));
figure(2)
mesh(X,Y,datas{5}.Fc(X,Y)+datas{6}.Fc(X,Y))
%%
for k=[1,5,6,8]
figure(k)
mesh(X,Y,datas{k}.Fc(X,Y)+datas{5}.Fc(X,Y)*0)
end

%mesh(X,Y,datas{3}.Fc(X,Y)./datas{4}.Fc(X,Y))
%% grid compearing
mesh(X,Y,datas{2}.Fc(X,Y)),set(gca,'ZDir','reverse')
mesh(X,Y,datas{2}.Fc(X,Y),smodels{2}.Ftc(X,Y)-datas{2}.Fc(X,Y)),set(gca,'ZDir','reverse'),colorbar
%mesh(X,Y,datas{2}.Fc(X,Y),smodels{2}.Ftc(X,Y)-smodels{1}.Ftc(X,Y)),set(gca,'ZDir','reverse'),colorbar
origo=[min(X(:)),min(Y(:))];
%axis([origo(1)+1.5e3,origo(1)+3e3,origo(2),origo(2)+6e3])

%% comparing all grids
figure(),clf,mesh(X,Y,Ftu(X,Y)),hold on,mesh(X,Y,abs(datas{2}.Fc(X,Y))),set(gca,'Zdir','reverse'),mesh(X,Y,smodels{1}.Ftc(X,Y))
%%
%plotting utsira grid from NPD with seismic data form Etor in ms
figure(),clf,mesh(X,Y,Ftu(X,Y)),hold on,mesh(X,Y,abs(datas{1}.Fc(X,Y))),set(gca,'Zdir','reverse')%,mesh(X,Y,smodels{2}.Ftc(X,Y))
%
%% comparing all grids
figure(),clf,mesh(X,Y,Ftu(X,Y).*0.92),hold on,mesh(X,Y,abs(datas{2}.Fc(X,Y))),set(gca,'Zdir','reverse'),mesh(X,Y,smodels{2}.Ftc(X,Y))


%% plot to be compeared with Ann Furre Sleipner 2014
aax=[4.3775e5 4.3975e5 6.47e6 6.474e6];
figure(1),clf,pcolor(X,Y,datas{7}.Fc(X,Y)),shading flat,caxis([0 11]),colorbar,axis equal,axis(aax)

%%
[ann ,map] = imread('data/Ann Furre_Sleipner 2014.tiff');
cmp=colormap();
%%
figure(),image(ann)
sz=size(ann);
mval=nan(sz(1),sz(2));
tol=0.25;
for i=1:sz(1)
    i
    for j=1:sz(2)
        val=squeeze(double(ann(i,j,:))/255)';
        d=sum(val);
        if((abs(d - 3) < tol) || abs(d)<tol)
            mval(i,j)=nan;
        else
            [v,ind]=min(sum(bsxfun(@minus,cmp,val).^2,2));
            mval(i,j)=ind;
        end
    end   
end
%%
cut=[40,925,97,520+10]
nmval=mval(cut(1):cut(2),cut(3):cut(4));
dx=10;dy1=10;dy2=20;
nmval(1:dx,:)=nan;
nmval(end-dx:end,:)=nan;
nmval(:,end-dy2:end)=nan;
nmval(:,1:dy1)=nan;
nmval=nmval*11/64;
x=linspace(aax(3),aax(4),cut(2)-cut(1)+1);
y=linspace(aax(1),aax(2),cut(4)-cut(3)+1);
[ya,xa]=meshgrid(x,y);
nmval=nmval';
%nmval=nmval(end:-1:1,:);
nmval=nmval(:,end:-1:1);
%xa=xa(end:-1:1,:);
figure(),pcolor(xa,ya,nmval);shading flat,colorbar,axis equal;axis tight
%%
mval=
pcolor(mval(43:920,100:519)),shading flat
%%
pcolor(mval),shading flat
