function y = interp2DTable(T, x, y)
    in_range_x = ( min(T.x)<=min(double(x)) && max(double(x)) <= max(T.x));
    in_range_y = ( min(T.y)<=min(double(y)) && max(double(y)) <= max(T.y));
    if(~(in_range_x && in_range_y))
       warning('value out of range of tensor table') 
    end
    
    if(isa(x,'ADI') || isa(y,'ADI'))
        xi=double(x);vi=double(y);
        [yi, dyidxi, dyidvi] = interp2DTableLoc(T, xi, vi);
        tmp={};
        if(isa(x,'ADI'))
         tmp{end+1} = lMultdiag(dyidxi, x.jac);
        end
        if(isa(y,'ADI'))
         tmp{end+1} = lMultdiag(dyidvi, y.jac);
        end
        if(numel(tmp)==1)
           y = ADI(yi, tmp{1});
           return
        elseif( numel(tmp)==2)   
            for i=1:numel(tmp{2})
                tmp{1}{i}=tmp{1}{i}+tmp{2}{i};
            end
        else
           error('Not possible')
        end                       
        y = ADI(yi, tmp{1});
        
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
    dyidxil = zeros(size(xi));
    dyidxir = zeros(size(xi));
    dwdvi   = zeros(size(xi));
end

%for k = 1:nreg
    xval=T.x;
    yval=T.y;
    %v = T{k}.key;  pos = T{k}.pos;
    v=yval;
    
    lims = v; lims(1)=-inf; lims(end)=inf;
    vi_k = xi;
    [bin,bin] = histc(vi, lims);                                     %#ok
    w = (vi-v(bin))./(v(bin+1)-v(bin));
    if compDer
        dwdvi = 1./(v(bin+1)-v(bin));
    end
    for tn = 1:numel(v)-1
        lns = tn:tn+1;
        tab = T.data(:,lns);
        %ixl = (bin==(tn-1));
        ixr = (bin==tn);
        if(any(ixr))
            yil(ixr) = fastInterpTable(xval, tab(:,1), vi_k(ixr));        
            yir(ixr) = fastInterpTable(xval, tab(:,2), vi_k(ixr));
        end
        if compDer
            dyidxil(ixr) = dinterpq1(xval, tab(:,1), vi_k(ixr));
            dyidxir(ixr) = dinterpq1(xval, tab(:,2), vi_k(ixr));
        end
    end
%end
yi = yil.*(1-w) + yir.*w;
if compDer
    dyidxi = dyidxil.*(1-w) + dyidxir.*w;
    dyidvi = -(yil-yir).*dwdvi;
end
end
