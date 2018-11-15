function x = getIcon(fn)
x = double(imread(fn))/255;
% set ~white to transparent
ix = sum(x,3) > 2.95;
ix = ix(:);
x([ix; ix; ix]) = nan;
end
