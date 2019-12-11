function [state]=FlowNMPFA(G,bc,fluid,W,OSflux,tol,maxiter)
[~,rho]=fluid.properties();
u=ones(G.cells.num,1)*mean([W(1).val,W(2).val]);
T=TransNMPFA(u);
[A,b]=AssemAb(T);
iter=0;
res=zeros(maxiter+1,1);
res(1)=norm(A*u-b);
while(res(iter+1)>tol*res(1)&&iter<maxiter)
    u=A\b;
    T=TransNMPFA(u);
    [A,b]=AssemAb(T);
    iter=iter+1;
    res(iter+1)=norm(A*u-b);
end
% u=A\b;T=TransNMPFA(u);
[flux,wellsol]=computeFlux(u,T);
state.pressure=u;state.flux=flux;state.wellSol=wellsol;
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
                
                r1=t1(3:end-1,2)'*(u(c1)-u(t1(3:end-1,1)));
                r1=r1+(t11-sum(t1(2:end-1,2)))*u(c1)-t1(end,2);
                r2=t2(3:end-1,2)'*(u(c2)-u(t2(3:end-1,1)));
                r2=r2+(t22-sum(t2(2:end-1,2)))*u(c2)-t2(end,2);
                
                eps=1e-6*max(abs([t1(:,end);t2(:,end)]));
                if((abs(r1)+abs(r2))<=eps)
                    mu1=0.5;mu2=0.5;
                else
                    mu1=abs(r2)/(abs(r1)+abs(r2));mu2=1-mu1;
                end
                
                if(r1*r2<=0)
                    T1=t1;T2=t2;
                    T1(1,2)=mu1*t12+mu2*t21+2*mu1*(t11-t12);
                    T1(2,2)=-(mu1*t12+mu2*t21);
                    T1(3:end,2)=-2*mu1.*T1(3:end,2);
                    T2(1,2)=mu1*t12+mu2*t21+2*mu2*(t22-t21);
                    T2(2,2)=-(mu1*t12+mu2*t21);
                    T2(3:end,2)=-2*mu2.*T2(3:end,2);
                else
                    trans=mu1*t12+mu2*t21;
                    T1=[c1 trans;c2 -trans;0 0];
                    T2=[c2 trans;c1 -trans;0 0];
                end
                T(i_face,1)={T1};
                T(i_face,2)={T2};
            else %boundary faces---------------
                ind=find(bc.face==i_face,1);
                if(strcmp(bc.type{ind},'pressure'))
                    t1=OSflux{i_face,1};
                    t1(2:end,2)=-t1(2:end,2);
                    T(i_face,1)={t1};
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
        
        
        for i=1:numel(W)
            if(strcmpi(W(i).type,'bhp'))
                pbh=W(i).val;dZ=W(i).dZ;
                for j=1:numel(W(i).cells)
                    mycell=W(i).cells(j);
                    I(k)=mycell;J(k)=mycell;V(k)=W(i).WI(j);k=k+1;
                    b(mycell)=b(mycell)+W(i).WI(j)*(pbh+rho*9.81*dZ(j));
                end
            else
                % write code here bababababababbababaabababababababababababababba
                error('code under development!')
            end
        end
        
        I(k:end)=[];J(k:end)=[];V(k:end)=[];
        A=sparse(I,J,V,nc,nc);
    end
%--------------------------------------------------------------------------
    function [flux,wellsol]=computeFlux(u,T)
        flux=zeros(G.faces.num,1);
        %----------------------------------------------
%         for i_face=1:G.faces.num
%             if(all(G.faces.neighbors(i_face,:)~=0))
%                 t1=T{i_face,1};t2=T{i_face,2};
%                 flux_left=t1(1:end-1,2)'*u(t1(1:end-1,1))+t1(end,2);
%                 flux_right=t2(1:end-1,2)'*u(t2(1:end-1,1))+t2(end,2);
%                 flux_left+flux_right
%                 flux(i_face)=0.5*(flux_left-flux_right);
%             else
%                 t1=T{i_face,1};
%                 flux(i_face)=t1(1:end-1,2)'*u(t1(1:end-1,1))+t1(end,2);
%             end
%         end
        %------------------------------------------------
        for i_face=1:G.faces.num
            t1=T{i_face,1};
            flux(i_face)=t1(1:end-1,2)'*u(t1(1:end-1,1))+t1(end,2);            
        end
        ind=G.faces.neighbors(:,1)==0;
        flux(ind)=-flux(ind);
              
        wellsol=repmat(struct('pressure',[],'flux',[]),[numel(W) 1]);
        for i=1:numel(W)
            if(strcmpi(W(i).type,'bhp'))
                pbh=W(i).val;dZ=W(i).dZ;
                wellsol(i).pressure=pbh+rho*9.81*dZ;
                wellsol(i).flux=W(i).WI.*(wellsol(i).pressure-u(W(i).cells));
            else
                error('code under development!');
                % write code here
            end
        end
    end
end