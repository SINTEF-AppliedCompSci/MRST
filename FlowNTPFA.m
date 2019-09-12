function [state]=FlowNTPFA(G,bc,fluid,W,OSflux,u0,tol,maxiter,varargin)

opt = struct('MatrixOutput', false, ...
             'src', []);
opt = merge_options(opt, varargin{:});

dispif(mrstVerbose, 'FlowNTPFA\n');
[mu,rho]=fluid.properties();
mu = mu(1);
rho = rho(1);

% Expand u0 if there are rate wells
numWR = getNoRateWells(W);
u0 = [u0; ones(numWR, 1)];
T=TransNTPFA(u0);
[A,b]=AssemAb(T,opt.src); 
iter=0;
res=zeros(maxiter+1,1);
bnorm=norm(b,inf);
res(1)=norm(A*u0-b,inf)/bnorm;
%res(1)=tol+1;
u=u0;
while(res(iter+1)>tol && iter<maxiter)
    dispif(mrstVerbose & mod(iter,10)==0, ['iter=',num2str(iter), ' res=', num2str(res(iter+1)), '\n'])
    u=A\b;
    T=TransNTPFA(u);
    [A,b]=AssemAb(T,opt.src);
    iter=iter+1;
    res(iter+1)=norm(A*u-b,inf)/bnorm;
end
dispif(mrstVerbose, ['iter=',num2str(iter), ' res=', num2str(res(iter+1)), '\n'])
if iter==maxiter
    warning(['Maximum number of iterations ', num2str(maxiter), ...
             ' reached. Residual=', num2str(res(end))]);
end

[flux,wellsol]=computeFlux(u,T);
state.pressure=u(1:G.cells.num);
state.flux=flux;
state.wellSol=wellsol;
state.iter=iter;
state.res=res(1:iter+1);

if opt.MatrixOutput
    state.A = A;
end
    
%--------------------------------------------------------------------------
function T=TransNTPFA(u)
    T=zeros(G.faces.num,2);
    for i_face=1:G.faces.num
        if(all(G.faces.neighbors(i_face,:)~=0))
            t1=OSflux{i_face,1};
            t2=OSflux{i_face,2};
            r1=t1(3:end-1,2)'*u(t1(3:end-1,1))+t1(end,2);
            r2=t2(3:end-1,2)'*u(t2(3:end-1,1))+t2(end,2);
            eps=1e-12*max(abs([t1(:,end);t2(:,end)]));
            if(abs(r1)<=eps),r1=0;end
            if(abs(r2)<=eps),r2=0;end
            
            if(abs(r1+r2)>eps)
                mu1=r2/(r1+r2);
                mu2=1-mu1;
            else
                mu1=0.5; 
                mu2=0.5;
            end
            T(i_face,1)=(mu1*t1(1,2)+mu2*t2(2,2))/mu;
            T(i_face,2)=(mu1*t1(2,2)+mu2*t2(1,2))/mu;
        else
            ind=find(bc.face==i_face,1);
            if(strcmpi(bc.type{ind},'pressure'))
                t1=OSflux{i_face,1};t2=OSflux{i_face,2};
                t11=t1(1,2);t12=t1(2,2);
                t22=t2(1,2);t21=t2(2,2);
                r1=t1(3:end-1,2)'*u(t1(3:end-1,1))+t1(end,2);
                r2=t2(end,2);
                eps=1e-12*max(abs([t1(:,end);t2(:,end)]));
                if(abs(r1)<=eps),r1=0;end
                if(abs(r2)<=eps),r2=0;end
                if(abs(r1+r2)>eps)
                    mu1=r2/(r1+r2);
                    mu2=1-mu1;
                else
                    mu1=0.5;
                    mu2=0.5;
                end
                T(i_face,1)=mu1*t11+mu2*t21; % Divide my mu?
                T(i_face,2)=(mu1*t12+mu2*t22)*bc.value{ind}(G.faces.centroids(i_face,:));
            else
                %T(i_face,2)=-G.faces.areas(i_face)*...
                %    bc.value{ind}(G.faces.centroids(i_face,:));
                T(i_face,2) = bc.value{ind}(G.faces.centroids(i_face,:));
            end
        end
    end
end
%--------------------------------------------------------------------------
function [A,b]=AssemAb(T, src)
    ncf=max(diff(G.cells.facePos));
    nc=G.cells.num;
    
    num_W_rate = getNoRateWells(W);
    nc = nc + num_W_rate;
    b=zeros(nc,1);
    [I,J,V]=deal(zeros(ncf*nc,1));k=1;
    for i_face=1:G.faces.num
        c1=G.faces.neighbors(i_face,1);
        c2=G.faces.neighbors(i_face,2);
        if(all([c1 c2]~=0))
            I(k)=c1;J(k)=c1;V(k)=T(i_face,1);k=k+1;
            I(k)=c1;J(k)=c2;V(k)=-T(i_face,2);k=k+1;
            I(k)=c2;J(k)=c2;V(k)=T(i_face,2);k=k+1;
            I(k)=c2;J(k)=c1;V(k)=-T(i_face,1);k=k+1;
        else
            c1=max(c1,c2);
            I(k)=c1;J(k)=c1;V(k)=T(i_face,1);k=k+1;
            b(c1)=b(c1)+T(i_face,2);
        end
    end
    
    %----------------------------------------------------------
    for i=1:numel(W)
        if(strcmpi(W(i).type,'bhp'))
            pbh=W(i).val;dZ=W(i).dZ;
            for j=1:numel(W(i).cells)
                mycell=W(i).cells(j);
                I(k)=mycell;J(k)=mycell;V(k)=W(i).WI(j)/mu;k=k+1;
                b(mycell)=b(mycell)+W(i).WI(j)/mu*(pbh+rho*9.81*dZ(j));
            end
        elseif(strcmpi(W(i).type,'rate'))
            rate=W(i).val;dZ=W(i).dZ;
            for j=1:numel(W(i).cells)
                mycell=W(i).cells(j);
                I(k)=mycell;J(k)=mycell;V(k)=W(i).WI(j)/mu;k=k+1;
                I(k)=mycell;J(k)=G.cells.num+i;V(k)=-W(i).WI(j)/mu;k=k+1;
                I(k)=G.cells.num+i;J(k)=G.cells.num+i;V(k)=W(i).WI(j)/mu;k=k+1;
                I(k)=G.cells.num+i;J(k)=mycell;V(k)=-W(i).WI(j)/mu;k=k+1;
                b(mycell)=b(mycell)+W(i).WI(j)/mu*(rho*9.81*dZ(j));
                b(G.cells.num+i)=b(G.cells.num+i)+W(i).WI(j)/mu*(rho*9.81*dZ(j));
            end
            b(G.cells.num+i)=b(G.cells.num+i)+rate;
        else
            % write code here bababababababbababaabababababababababababababba
            error('code under development!')
        end
    end
    %-------------------------------------------------------
    I(k:end)=[];J(k:end)=[];V(k:end)=[];
    A=sparse(I,J,V,nc,nc);
    if ~isempty(src)
        b(src.cell(:))=b(src.cell(:))+src.rate(:);
    end
end
%--------------------------------------------------------------------------
function [flux,wellsol]=computeFlux(u,T)
    flux=zeros(G.faces.num,1);
    ind=all(G.faces.neighbors~=0,2);
    c1=G.faces.neighbors(ind,1);c2=G.faces.neighbors(ind,2);
    flux(ind)=T(ind,1).*u(c1)-T(ind,2).*u(c2);
    c1=max(G.faces.neighbors(~ind,:),[],2);
    flux(~ind)=T(~ind,1).*u(c1)-T(~ind,2);
    ind=G.faces.neighbors(:,1)==0;
    flux(ind)=-flux(ind);
    
    wellsol=repmat(struct('pressure',[],'flux',[]),[numel(W) 1]);
    for i=1:numel(W)
        if(strcmpi(W(i).type,'bhp'))
            pbh=W(i).val;
            dZ=W(i).dZ;
            wellsol(i).pressure=pbh+rho*9.81*dZ;
            wellsol(i).flux=W(i).WI./mu.*(wellsol(i).pressure-u(W(i).cells));
        elseif(strcmpi(W(i).type,'rate'))
            rate=W(i).val;
            %dZ=W(i).dZ;
            wellsol(i).pressure=u(G.cells.num+i);
            wellsol(i).flux=W(i).WI./mu.*(wellsol(i).pressure-u(W(i).cells)+rate);
        else
            error('code under development!');
            % write code here babbabababaababababababababababababababababababab
        end
    end
end
end 

function num_W_rate = getNoRateWells(W)
    num_W_rate = 0;
    for i=1:numel(W)
        if(strcmpi(W(i).type,'rate'))
            num_W_rate = num_W_rate+1;
        end
    end
end