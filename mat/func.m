function f = func(X,x1,y1)

f = 0;
for i = 1:numel(x1)
    f = f + 1/sqrt(2*delta*pi).*exp(-((X(:,1)-x1').^2 + (X(:,2)-y1').^2)./(2*delta));
end

end
