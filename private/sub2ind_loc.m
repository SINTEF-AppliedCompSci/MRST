function ndx = sub2ind_loc(siz,varargin)


%{ 
Copyright 2014 Jan Nordbotten jan.nordbotten@math.uib.no
%} 
siz = double(siz); 

%Compute linear indices
k = [1 cumprod(siz(1:end-1))];
ndx = 1;
s = size(varargin{1}); %For size comparison
for i = 1:length(siz),
    v = varargin{i};
    ndx = ndx + (v-1)*k(i);
end
