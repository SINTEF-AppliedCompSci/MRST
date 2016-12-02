function [uu, extra] = LMIM_linElast(G, C, el_bc, load, varargin)
opt = struct('solvetype','direct');
opts =  merge_options(opt,varargin{:});
opt=struct('stabterm','inter','as_type','CN_type','invertBlocks', 'mex');
tic;
lsmim = localStressMimetic(G, C,'theta',0.3,'invertBlocks',opt.invertBlocks,'stabterm',opt.stabterm,'as_type',opt.as_type);
vec=[1 0];
%profile viewer

mm='direct';
switch mm
    case 'mixed'
        ngam=size(lsmim.mixed.gama1,1);
        rhs=[zeros(size(lsmim.mixed.B,1),1);...
            reshape(bsxfun(@times,load(G.cells.centroids),G.cells.volumes)',[],1);...
            zeros(ngam,1)];      
        mm1=[lsmim.mixed.CT',lsmim.mixed.gama1'];
        mm2=[lsmim.mixed.CT',lsmim.mixed.gama2'];
        A=[lsmim.mixed.B,mm1;mm2',sparse(size(mm1,2),size(mm1,2))];
        u=A\rhs;
        utmp=u(size(lsmim.mixed.B,1)+1:size(lsmim.mixed.B,1)+G.cells.num*G.griddim);
        uu=-reshape(utmp,2,[])';
        if(ngam==G.nodes.num)
            figure(),plotNodeData(G,u(end-ngam:end)),colorbar
        end
    case 'direct'
        u=lsmim.direct.A\reshape(bsxfun(@times,load(G.cells.centroids),G.cells.volumes)',[],1);
        uu=reshape(u,2,[])';
    otherwise 
        error()
end

extra=lsmim;
end