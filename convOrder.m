function [r] = convOrder(nc,er,dim)
%Compute convergence order of pressure/flux solutions

r=zeros(size(er));
for i=2:numel(nc)
    for j=1:size(er,2)
        r(i,j)=-dim*log(er(i,j)/er(i-1,j))/log(nc(i)/nc(i-1));
    end
end

