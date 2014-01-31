function info=refineInfo(info_old,ref2);
%ref=2;
%ref=2^ref2;
ref=ref2;
info=info_old;
meta_old=info_old.meta;

% refine meta
meta=meta_old;
meta.ncols=(meta_old.ncols-1)*ref+1;
meta.nrows=(meta_old.nrows-1)*ref+1;
meta.cellsize=meta_old.cellsize/ref;
meta.dims=[meta_old.dims-1]*ref+1;

data_old=info_old.data;


%data = interp2(data_old,ref2);
%[Xn,Yn]=meshgrid(linspace(1,meta_old,nrows,meta.nrows),linspace(1,meta_old,ncols,meta.ncols),
[Xn,Yn]=meshgrid(linspace(1,meta_old.ncols,meta.ncols),linspace(1,meta_old.nrows,meta.nrows));%linspace(1,meta_old,ncols,meta.ncols),
data = interp2(data_old,Xn,Yn);
clear Xn;
clear Yn;
info.meta=meta;
info.data=data;
end