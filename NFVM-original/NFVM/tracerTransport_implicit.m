function [state] = tracerTransport_implicit(G,rock,W,state,dt,nstep)
%Solve the passive tracer transport equation implicity
%   state.cwell - contrentation at the prodution well
%   state.cres - distribution of concentration throught the reservoir

nc=G.cells.num;
nf=G.faces.num;
nw=numel(W);
ncf=100;
state.cwell=zeros(nstep,1);
state.cres=zeros(nc,nstep+1);

pv=G.cells.volumes.*rock.poro;
for t=1:nstep
    b=zeros(nc,1);
    [I,J,V]=deal(zeros(ncf*nc,1));k=1;
    
    for i=1:nc
        I(k)=i;J(k)=i;V(k)=pv(i)/dt;k=k+1;
        b(i)=b(i)+pv(i)/dt*state.cres(i,t);
    end
    
    for i=1:nf
        c1=G.faces.neighbors(i,1);
        c2=G.faces.neighbors(i,2);
        if(all([c1 c2]~=0))
            if(state.flux(i)>0)
                I(k)=c1;J(k)=c1;V(k)=state.flux(i);k=k+1;
                I(k)=c2;J(k)=c1;V(k)=-state.flux(i);k=k+1;
            else
                I(k)=c1;J(k)=c2;V(k)=state.flux(i);k=k+1;
                I(k)=c2;J(k)=c2;V(k)=-state.flux(i);k=k+1;
            end
        end
    end
    
    for i=1:nw
        if(strcmpi(W(i).name,'I'))
            b(W(i).cells)=b(W(i).cells)+state.wellSol(i).flux;
        else
            for j=1:numel(W(i).cells)
                mycell=W(i).cells(j);
                I(k)=mycell;J(k)=mycell;V(k)=-state.wellSol(i).flux(j);k=k+1;
            end
        end
    end
    I(k:end)=[];J(k:end)=[];V(k:end)=[];
    A=sparse(I,J,V,nc,nc);
    x=A\b;
    state.cres(:,t+1)=x;
    
    for i=1:nw
        if(strcmpi(W(i).name,'P'))
            mycells=W(i).cells;
            myconcentration=x(mycells);
            myflux=-state.wellSol(i).flux;
            state.cwell(t)=dot(myconcentration,myflux)/sum(myflux);
            break;
        end
    end
    
end
end

