function h=cplot3(x,y,z,c)
x = reshape(x,[],1);
y = reshape(y,[],1);
z = reshape(z,[],1);
c = reshape(c,[],1);
h = surface([x x],[y y],[z z],[c c],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
end

