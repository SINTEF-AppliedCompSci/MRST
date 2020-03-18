function state=PicardNMPFA(G,bc,src,OSflux,u0,tol,maxiter)
% Picard iteration for nonlinear MPFA method
T=TransNMPFA(u0);
[A,b]=AssemAb(T);
iter=0;
res=zeros(maxiter+1,1);
res(1)=norm(A*u0-b,inf);
while(res(iter+1)>tol*res(1)&&iter<maxiter)
    u=A\b;
%     uu=bicgstab(A,b,1e-12,100);
%     uuu=gmres(A,b,[],1e-12,100);
    T=TransNMPFA(u);
    [A,b]=AssemAb(T);
    iter=iter+1;
    res(iter+1)=norm(A*u-b,inf);
end
flux=computeFlux(u,T);
state.pressure=u;state.flux=flux;
state.iter=iter;state.res=res(1:iter+1);
%--------------------------------------------------------------------------
    function T=TransNMPFA(u)
        T=cell(G.faces.num,2);
        for i_face=1:G.faces.num
            if(all(G.faces.neighbors(i_face,:)~=0)) % internal face
                c1=G.faces.neighbors(i_face,1);
                c2=G.faces.neighbors(i_face,2);
                t1=OSflux{i_face,1};
                t2=OSflux{i_face,2};
                t11=t1(1,2);t12=t1(2,2);
                t22=t2(1,2);t21=t2(2,2);
                r1=t1(3:end-1,2)'*(u(c1)-u(t1(3:end-1,1)))+...
                    (t11-sum(t1(2:end-1,2)))*u(c1)-t1(end,2);
                r2=t2(3:end-1,2)'*(u(c2)-u(t2(3:end-1,1)))+...
                    (t22-sum(t2(2:end-1,2)))*u(c2)-t2(end,2);               
                eps=1e-12*max(abs([t1(:,end);t2(:,end)]));
                if(abs(r1)<=eps),r1=0;end
                if(abs(r2)<=eps),r2=0;end
                %----------------------------------------------------------
                if((abs(r1)+abs(r2))<=eps)
                    mu1=0.5;mu2=0.5;
                else
                    mu1=abs(r2)/(abs(r1)+abs(r2));mu2=1-mu1;
                end
                trans=mu1*t12+mu2*t21;
                if(r1*r2<0)
                    T1=t1;T2=t2;
                    T1(1,2)=trans+2*mu1*(t11-t12);
                    T1(2,2)=-trans;
                    T1(3:end,2)=-2*mu1.*T1(3:end,2);
                    T2(1,2)=trans+2*mu2*(t22-t21);
                    T2(2,2)=-trans;
                    T2(3:end,2)=-2*mu2.*T2(3:end,2);
                else
                    T1=[c1 trans;c2 -trans;0 0];
                    T2=[c2 trans;c1 -trans;0 0];
                end
                %----------------------------------------------------------
%                 if(r1==0&&r2==0)
%                     mu1=0.5;mu2=0.5;beta1=1;beta2=1;
%                     trans=mu1*t12+mu2*t21;
%                 elseif(r1==0)
%                     mu1=1;mu2=0;beta1=1;beta2=0;
%                     trans=0.5*t12+0.5*t21;
%                 elseif(r2==0)
%                     mu1=0;mu2=1;beta1=0;beta2=1;
%                     trans=0.5*t12+0.5*t21;
%                 else
%                     mu1=abs(r2)/(abs(r1)+abs(r2));mu2=1-mu1;
%                     beta1=mu1*(1-sign(r1*r2));
%                     beta2=mu2*(1-sign(r1*r2));
%                 end
%                 T1=t1;T2=t2;
%                 T1(1,2)=trans+beta1*(t11-t12);
%                 T1(2,2)=-trans;
%                 T1(3:end,2)=-beta1*T1(3:end,2);
%                 T2(1,2)=trans+beta2*(t22-t21);
%                 T2(2,2)=-trans;
%                 T2(3:end,2)=-beta2*T2(3:end,2);
                %----------------------------------------------------------
%                 mu1=(abs(r2)+eps)/(abs(r1)+abs(r2)+2*eps);
%                 mu2=1-mu1;
%                 beta1=mu1*(1-sign(r1*r2));
%                 beta2=mu2*(1-sign(r1*r2));
%                 trans=mu1*t12+mu2*t21;
%                 T1=t1;T2=t2;
%                 T1(1,2)=trans+beta1*(t11-t12);
%                 T1(2,2)=-trans;
%                 T1(3:end,2)=-beta1*T1(3:end,2);
%                 T2(1,2)=trans+beta2*(t22-t21);
%                 T2(2,2)=-trans;
%                 T2(3:end,2)=-beta2*T2(3:end,2);
                %----------------------------------------------------------
                T(i_face,1)={T1};
                T(i_face,2)={T2};
            else %boundary faces---------------
                ind=find(bc.face==i_face,1);
                if(strcmp(bc.type{ind},'pressure'))
                    t1=OSflux{i_face,1};t2=OSflux{i_face,2};
                    t11=t1(1,2);t12=t1(2,2);
                    t22=t2(1,2);t21=t2(2,2);
                    r1=t1(3:end-1,2)'*u(t1(3:end-1,1))+t1(end,2);
                    r2=t2(end,2);
                    eps=1e-12*max(abs([t1(:,end);t2(:,end)]));
                    if(abs(r1)<=eps),r1=0;end
                    if(abs(r2)<=eps),r2=0;end
                    if(abs(r1+r2)>eps)
                        mu1=r2/(r1+r2);mu2=1-mu1;
                    else
                        mu1=0.5;mu2=0.5;
                    end
                    c1=max(G.faces.neighbors(i_face,:));
                    trans=mu1*t11+mu2*t21;
                    const=(mu1*t12+mu2*t22)*bc.value{ind}(G.faces.centroids(i_face,:));
                    T(i_face,1)={[c1 trans;0 -const]};
                else
                    gN=G.faces.areas(i_face)*...
                        bc.value{ind}(G.faces.centroids(i_face,:));
                    T(i_face,1)={[0 gN]};
                end
            end
        end
    end
%--------------------------------------------------------------------------
    function [A,b]=AssemAb(T)
        ncf=50;
        nc=G.cells.num;
        b=zeros(nc,1);
        [I,J,V]=deal(zeros(ncf*nc,1));k=1;
        for i_face=1:G.faces.num
            c1=G.faces.neighbors(i_face,1);
            c2=G.faces.neighbors(i_face,2);
            if(all([c1 c2]~=0))
                t1=T{i_face,1};t2=T{i_face,2};
                I(k:k+size(t1,1)-2)=c1;
                J(k:k+size(t1,1)-2)=t1(1:end-1,1);
                V(k:k+size(t1,1)-2)=t1(1:end-1,2);
                k=k+size(t1,1)-1;
                b(c1)=b(c1)-t1(end,2);
                I(k:k+size(t2,1)-2)=c2;
                J(k:k+size(t2,1)-2)=t2(1:end-1,1);
                V(k:k+size(t2,1)-2)=t2(1:end-1,2);
                k=k+size(t2,1)-1;
                b(c2)=b(c2)-t2(end,2);
            else
                c1=max(c1,c2);
                t1=T{i_face,1};
                I(k:k+size(t1,1)-2)=c1;
                J(k:k+size(t1,1)-2)=t1(1:end-1,1);
                V(k:k+size(t1,1)-2)=t1(1:end-1,2);
                k=k+size(t1,1)-1;
                b(c1)=b(c1)-t1(end,2);
            end
        end
        I(k:end)=[];J(k:end)=[];V(k:end)=[];
        A=sparse(I,J,V,nc,nc);
        b(src.cell(:))=b(src.cell(:))+src.rate(:);
    end
%--------------------------------------------------------------------------
    function flux=computeFlux(u,T)
        flux=zeros(G.faces.num,1);
        for i_face=1:G.faces.num
            t1=T{i_face,1};
            flux(i_face)=t1(1:end-1,2)'*u(t1(1:end-1,1))+t1(end,2);
        end
        ind=G.faces.neighbors(:,1)==0;
        flux(ind)=-flux(ind);
    end
end