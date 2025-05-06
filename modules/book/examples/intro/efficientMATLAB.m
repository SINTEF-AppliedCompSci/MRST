%%
n = 5000000;
pt = randn(n,3);

%%
tic
I = sum(bsxfun(@times, pt>0, [1 2 4]),2)+1;
num = accumarray(I,1);
toc

%%
tic
avg = zeros(8,3);
for j=1:size(pt,1)
    quad = sum((pt(j,:)>0).*[1 2 4])+1;
    avg(quad,:) = avg(quad,:)+pt(j,:);
end
avg = bsxfun(@rdivide, avg, num);
toc

%%
tic
avg = bsxfun(@rdivide, sparse(I,1:n,1)*pt, accumarray(I,1));
toc

%%
tic
avg = sparse(I,1:n,1)*[pt, ones(n,1)];
avg = bsxfun(@rdivide, avg(:,1:end-1), avg(:,end));
toc