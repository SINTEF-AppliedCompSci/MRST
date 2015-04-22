%%
ff=exp(-((xx-xc(400))/(0.3)).^2);ff=ff/sum(ff);plot(ff)
%%
xc=Gt.cells.centroids(:,1)/1e3;
z=Gt.cells.z;
z_org=Gt_org.cells.z;
z_new=max(z_org)*ones(size(z_org));
%assert(z_org(1)==z_new(1));
%%
for i=1:numel(z_new)-1
    z_new(i)=max(z_org(i:end));
end
z_new(end)=max(z_new(end-1),z_org(end))
plot(xc,[z,z_org,z_new])
%plot(xc,[z_org,z_new])
legend('org','new')
shirt=100
pp=numel(xc)-700;
axis([xc(pp-shirt) xc(pp) min(z_org(pp-shirt:pp)) max(z_org(pp-shirt:pp))])
plot(xc,filter2(ff,z_new-z_org))
h_trap=filter2(ff,z_new-z_org)