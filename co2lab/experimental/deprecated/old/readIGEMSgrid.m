function [Gt,rock2D,SVE,rock,invert_axis,G]=readIGEMSgrid(mygrid,mymod)
%   [Gt,rock2D,SVE,rock,invert_axis,G]=readIGEMSgrid(mygrid,mymod)
%   function which read and process igems eclipse grids or read
%   preprocessed  values.
%   mygrid - name of grid file mygrid.GRDECL
%   mymod  - modification mode, see code for definition 
%            > edit readIGEMSgrid.m
   igems_dir = mfilename('fullpath');
   mind=find(igems_dir=='/');
   igems_dir=igems_dir(1:mind(end));
   run_name=[mygrid,'_',mymod];
   try
      %error()
      %aa=load(['/data/igems/30times60km/matlab/',run_name,'.mat']);
      aa=load(fullfile(igems_dir,'data/matlab/',[run_name,'.mat']));
      Gt=aa.Gt;
      %WVE=aa.WVE;
      rock2D=aa.rock2D;
      %rate=aa.rate;
      %W=aa.W;
      SVE=aa.SVE;
      rock=aa.rock;
      invert_axis=aa.invert_axis;
      if(nargout==6)
         G=aa.G;
         assert(size(rock.perm,1)==G.cells.num);  
      end
   catch 
      disp(['Reading and processing grid for ',mygrid,'  ',mymod])
      disp(' -> Reading data'); tic
      grdecl_filename=fullfile(igems_dir,'data/eclipsegrids/',[mygrid,'.GRDECL']);
      grdecl = readGRDECL(grdecl_filename); toc
      disp(' -> Processing grid'); tic
      switch mymod
         case 'original'
            disp('Use original grid')
         case 'z1'
            %grdecl_old=grdecl;
            grdecl.PERMY=grdecl.PERMX;
            grdecl.PERMZ=grdecl.PERMX;
            grdecl=coarseGrdecl(grdecl,[1 1 grdecl.cartDims(3)]);
            grdecl.ACTNUM=int32(grdecl.ACTNUM);
        case 'z1_uniform'
            %grdecl_old=grdecl;
            grdecl=coarseGrdecl(grdecl,[1 1 grdecl.cartDims(3)],'only_grid',true);
            grdecl.PERMX=ones(prod(grdecl.cartDims),1)*500;%*milli*darcy;
            grdecl.PERMY=grdecl.PERMX;
            grdecl.PERMZ=grdecl.PERMX;
            grdecl.PORO=ones(prod(grdecl.cartDims),1)*0.2;
            grdecl.ACTNUM=int32(grdecl.ACTNUM);
         case 'z1_coarse22'
            %grdecl_old=grdecl;
            grdecl=coarseGrdecl(grdecl,[2 2 grdecl.cartDims(3)]);
            grdecl.PERMY=grdecl.PERMX;
            grdecl.PERMZ=grdecl.PERMX;
            grdecl.ACTNUM=int32(grdecl.ACTNUM);           
         case 'z1_coarse22_uniform'
            %grdecl_old=grdecl;
            grdecl=coarseGrdecl(grdecl,[2 2 grdecl.cartDims(3)],'only_grid',true);
            grdecl.PERMX=ones(prod(grdecl.cartDims),1)*500;%*milli*darcy;
            grdecl.PERMY=grdecl.PERMX;
            grdecl.PERMZ=grdecl.PERMX;
            grdecl.PORO=ones(prod(grdecl.cartDims),1)*0.2;
            grdecl.ACTNUM=int32(grdecl.ACTNUM);
         case 'z1_coarse44'
            %grdecl_old=grdecl;
            grdecl.PERMY=grdecl.PERMX;
            grdecl.PERMZ=grdecl.PERMX;
            grdecl=coarseGrdecl(grdecl,[4 4 grdecl.cartDims(3)]);
            grdecl.ACTNUM=int32(grdecl.ACTNUM);
         case 'coarse44'
            %grdecl_old=grdecl;
            grdecl.PERMY=grdecl.PERMX;
            grdecl.PERMZ=grdecl.PERMX;
            grdecl=coarseGrdecl(grdecl,[4 4 1]);
            grdecl.ACTNUM=int32(grdecl.ACTNUM);
         case 'z1_coarse44_uniform'
            %grdecl_old=grdecl;
            grdecl=coarseGrdecl(grdecl,[4 4 grdecl.cartDims(3)],'only_grid',true);
            grdecl.PERMX=ones(prod(grdecl.cartDims),1)*500;%*milli*darcy;
            grdecl.PERMY=grdecl.PERMX;
            grdecl.PERMZ=grdecl.PERMX;
            grdecl.PORO=ones(prod(grdecl.cartDims),1)*0.2;
            grdecl.ACTNUM=int32(grdecl.ACTNUM);
         case 'z1_ref22_uniform'
            %grdecl_old=grdecl;
            %grdecl_org=grdecl;
            grdecl=coarseGrdecl(grdecl,[1 1 grdecl.cartDims(3)],'only_grid',true);
            grdecl=refineGrdecl(grdecl,[2 2 1]);
            grdecl.PERMX=ones(prod(grdecl.cartDims),1)*500;%*milli*darcy;
            grdecl.PERMY=grdecl.PERMX;
            grdecl.PERMZ=grdecl.PERMX;
            grdecl.PORO=ones(prod(grdecl.cartDims),1)*0.2;
            grdecl.ACTNUM=int32(grdecl.ACTNUM);
         otherwise            
            error(['no such mymod: ', mymod])
      end
      %grdecl.COORD(2:3:end)=-grdecl.COORD(2:3:end);
      [X,Y,Z]  = buildCornerPtNodes(grdecl);
      V1 = [X(end,1,1) - X(1,1,1), Y(end,1,1) - Y(1,1,1), 0                  ];
      V2 = [X(1,end,1) - X(1,1,1), Y(1,end,1) - Y(1,1,1), 0                  ];
      V3 = [0,                     0,                     Z(1,1,end)-Z(1,1,1)];
      invert_axis=false;
      if(cross(V1,V2)*V3' < 0)
         grdecl.COORD(2:3:end)=-grdecl.COORD(2:3:end);
         invert_axis=true;
      end
      G=processgrid(grdecl);
      %G=processGRDECL(grdecl);
      rock = grdecl2Rock(grdecl,G.cells.indexMap);
      %rock = initEclipseRock(grdecl,)
      %clear grdecl;
      tic;G = computeGeometry(G); toc
      rock.perm = convertFrom(rock.perm, milli*darcy);
      assert(size(rock.perm,1)==G.cells.num);
      % Top-surface grid
      disp(' -> Constructing top-surface grid'); tic
      [Gt, G] = topSurfaceGrid(G); toc
      rock.perm=rock.perm(:,1);
      rock2D  = averageRock(rock, Gt);
      %initilize innerproduct
      disp(' -> Initialising solvers'); tic
      %SVE = computeMimeticIPVE(Gt, rock2D, 'Innerproduct','ip_simple'); toc
      disp(' -> Not Initialising solvers'); tic
      SVE = [];
      if(~isdir(fullfile(igems_dir,'data/matlab')))
         mkdir(fullfile(igems_dir,'data/matlab'));
      end
      save(fullfile(igems_dir,'data/matlab/',[run_name,'.mat']),'-v7.3')
   end
