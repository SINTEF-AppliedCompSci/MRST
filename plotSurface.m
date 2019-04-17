function [ ] = plotSurface(coord,u,N,opt,ititle )
%Plot contour and surf plot of pressure solution
%   coord: coordinate of data points
%   u: pressure value
%   N: number of grid points
%   opt: 'griddata' - use MATLAB griddata function
%        'trisurf' - use MATLAB trisurf function
%        'Direct' - use user-provided grid

switch opt
    case 'griddata'
        x=coord(:,1);y=coord(:,2);
        [xq,yq] = meshgrid(linspace(min(x),max(x),N(1)),...
            linspace(min(y),max(y),N(2)));
        vq=griddata(x,y,u,xq,yq,'v4');
        figure,contourf(xq,yq,vq,20)
%         colorbar('southoutside')
        colorbar
        caxis([min(u) max(u)]);
        axis image off;title(ititle)
        figure,surf(xq,yq,vq)
        colorbar,caxis([min(u) max(u)]);title(ititle)
    case 'direct'
        assert(numel(u)==prod(N));
        x=reshape(coord(:,1),N(1),N(2));x=x';
        y=reshape(coord(:,2),N(1),N(2));y=y';
        p=reshape(u,N(1),N(2));p=p';
        figure,contourf(x,y,p,20);
%         colorbar('southoutside')
        colorbar;
        caxis([min(u) max(u)]);
        axis image off,title(ititle)
        figure,surf(x,y,p);
        colorbar,caxis([min(u) max(u)]);title(ititle);
    case 'trisurf'
        x=coord(:,1);y=coord(:,2);
        [xq,yq] = meshgrid(linspace(min(x),max(x),N(1)),...
            linspace(min(y),max(y),N(2)));
        vq=griddata(x,y,u,xq,yq,'v4');
        figure,contourf(xq,yq,vq,20)
        colorbar('southoutside'),caxis([min(u) max(u)]);
        axis image off;title(ititle)
        t=delaunay(x,y);
        figure,trisurf(t,x,y,u);title(ititle);
        colorbar;caxis([min(u) max(u)]);
    otherwise
        error('Unknown plot method!');
end
end

