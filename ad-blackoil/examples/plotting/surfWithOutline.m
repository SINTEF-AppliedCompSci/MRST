function surfWithOutline(x, y, z)
%SURFWITHOUTLINE - standard surf plot with black outline of perimeter

surf(x,y,z)
flag = false;
if ~ishold, flag=true; hold on; end
plot3(x(:,[1 end]), y(:,[1 end]), z(:,[1 end]),'-k','LineWidth',2)
plot3(x([1 end],:)', y([1 end],:)', z([1 end],:)','-k','LineWidth',2)
if flag, hold off, end
end

