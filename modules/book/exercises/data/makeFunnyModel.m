%% Make funny data set
mrstModule add libgeometry;
img = imread('mickey-small.png'); img=double(img(:,:,1));
img = img/max(img(:));
dim = [size(img) 1] +[26 15 9];
G = mcomputeGeometry(cartGrid(dim));
I = ones(dim);
I(14:end-13,10:end-6,end-3:end-3) = img; I = I(:);
K = logNormLayers(dim, repmat(400,dim(3),1)); disp([max(K(:)) min(K(:))])
K(I~=1) = 250+150*I(I~=1);