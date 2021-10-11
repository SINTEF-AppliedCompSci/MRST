function [sigm,evec]=calStressEigsVEM(G,stress,varargin)
% calculate eigen values of stress tensor and the basis
% 
% transform to traditional voits notation (for stress this has no factors)

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
%}

opt=struct('use_c',false);
opt = merge_options(opt,varargin{:}); 
if(G.griddim==2)
   stresstens=[stress(:,1),stress(:,3),...
               stress(:,3),stress(:,2)];
else
   %stresstens=[stress(:,1),stress(:,end),stress(:,end-1),...
   %            stress(:,end),stress(:,2),stress(:,end-2),...
   %            stress(:,end-1),stress(:,end-2),stress(:,3)];
   stresstens=[stress(:,1),stress(:,end),stress(:,end-1),...
               stress(:,end),stress(:,2),stress(:,end-2),...
               stress(:,end-1),stress(:,end-2),stress(:,3)];        
   warning('check consistency of stress tensor');
end
if(opt.use_c)
    [eigc,vec]=multiSymmEig(reshape(stresstens,[],1));
    sigm=reshape(eigc,G.griddim,[])';
    evec=reshape(vec,G.griddim.^2,[])';
else
    sigm=nan(G.cells.num,G.griddim);
    evec=nan(G.cells.num,G.griddim.^2);
    for i=1:G.cells.num
        lstress=reshape(stresstens(i,:),G.griddim,G.griddim);
        [le,lsigm]=eig(lstress);
        lsigm=diag(lsigm);
      %{
      % NEVER EVER use exept for plotting  
      if(abs(le(1,1))<abs(le(1,2)))
         le(:,[1,2])=le(:,[2,1]);
         lsigm([1:2])=lsigm([2,1]);
      end
      % hack
      if(le(1,1)<0)
          le=-le;
      end
      if(det(le)<0)
          le(:,2)=-le(:,2);
      end
        %}
        
        sigm(i,:)=lsigm;
        evec(i,:)=le(:);
    end
end
end
