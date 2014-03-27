function y = interp2DTable(T, x, y)
    in_range_x = ( min(T.x)<=min(x) && max(x) <= max(T.x));
    in_range_y = ( min(T.y)<=min(y) && max(y) <= max(T.y));
    if(~(in_range_x && in_range_y))
       warning('value out of range of tensor table') 
    end
    
    if(isa(x,'ADI') || isa(y,'ADI'))
        xi=double(x);vi=double(y);
        [yi, dyidxi, dyidvi] = interp2DTableLoc(T, xi, vi);
         tmp1 = lMultdiag(dyidxi, x.jac);
         tmp2 = lMultdiag(dyidvi, y.jac);
            for i=1:numel(tmp2)
                tmp1{i}=tmp1{i}+tmp2{i};
            end                       
            y = ADI(yi, tmp1);
        
    else
       y = interp2DTableLoc(T, x, y);
    end
end

function [yi, dyidxi, dyidvi] = interp2DTableLoc(T, xi, vi)
compDer = (nargout>1);

%w   = zeros(size(xi));
yil = zeros(size(xi));
yir = zeros(size(xi));
if compDer
    dyidvil = zeros(size(xi));
    dyidvir = zeros(size(xi));
    dwdvi   = zeros(size(xi));
end

%for k = 1:nreg
    xval=T.x;
    yval=T.y;
    %v = T{k}.key;  pos = T{k}.pos;
    v=xval;
    
    lims = v; lims(1)=-inf; lims(end)=inf;
    vi_k = vi;
    [bin,bin] = histc(xi, lims);                                     %#ok
    w = (xi-v(bin))./(v(bin+1)-v(bin));
    if compDer
        dwdvi = 1./(v(bin+1)-v(bin));
    end
    for tn = 1:numel(v)-1
        lns = tn:tn+1;
        tab = T.data(:,lns);
        %ixl = (bin==(tn-1));
        ixr = (bin==tn);
        if(any(ixr))
            yil(ixr) = interpTable(yval, tab(:,1), vi_k(ixr));        
            yir(ixr) = interpTable(yval, tab(:,2), vi_k(ixr));
        end
        if compDer
            dyidvil(ixr) = dinterpTable(yval, tab(:,1), vi(ixr));
            dyidvir(ixr) = dinterpTable(yval, tab(:,2), vi(ixr));
        end
    end
%end
yi = yil.*(1-w) + yir.*w;
if compDer
    dyidvi = dyidvil.*(1-w) + dyidvir.*w;
    dyidxi = -(yil-yir).*dwdvi;
end
end
