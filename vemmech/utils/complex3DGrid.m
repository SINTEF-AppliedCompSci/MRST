function [G,G_org]=complex3DGrid(opt,mycase)
G_org=[];
switch mycase
    case 'box'        
        G=squareGrid(opt.cartDims,opt.L,'grid_type','cartgrid','disturb',opt.disturb);
        G_old=computeGeometry(G);
        %face=G.faces.num-13 %3 3 1 mitt surface
        if(opt.triangulate)
            %face=1:G.faces.num;
            face=unique(G.cells.faces(any(bsxfun(@eq,G.cells.faces(:,2),[5,6]),2),1));
            G=triangulateFaces(G,[face]');
        end
        G=computeGeometry(G);
        G=mrstGridWithFullMappings(G)
    case 'grdecl'        
        %G=squareGrid(opt.cartDims,opt.L,'grid_type','cartgrid','disturb',opt.disturb);
        grdecl = simpleGrdecl([2, 1, 2]*ceil((1e3).^(1/griddim)), 0.15);
        G = processGRDECL(grdecl);
        G=mrstGridWithFullMappings(G);
        G=computeGeometry(G);
        %%
    case 'sbed'
        % extream 3D test case
        % grdecl=readGRDECL('C:\Users\hnil\Documents\BITBUCKET\mrst-mechanics-fractures\test_cases\data\WaveIrregular.grdecl');
        grdecl=readGRDECL('/home/xavier/Matlab/Projects/project-mechanics-fractures/test_cases/data/WaveIrregular.grdecl');
        grdecl.ZCORN=min(grdecl.ZCORN,5)-1;
        grdecl.ZCORN=max(grdecl.ZCORN,1)-1;
        grdecl_c=cutGrdecl(grdecl,[1 15;1 15; 1 333]);
        G=processGRDECL(grdecl_c,'Tolerance',opt.gtol);
        if(opt.triangulate)
            faces=unique(G.cells.faces(any(bsxfun(@eq,G.cells.faces(:,2),[5,6]),2),1));
            %faces=1:G.faces.num;
            %
         G=triangulateFaces(G,faces');
        end
        figure(),plotGrid(G)
        %%
    case 'model3'
        % extream 3D test case
        %%       
        grdecl=readGRDECL('/home/hnil/heim/SVN/simsvn/data/original_datafiles/test_cases/model3.grdecl')
        G=processGRDECL(grdecl);        
        G=mrstGridWithFullMappings(G);
        G=computeGeometry(G);
        figure()
        clf,plotGrid(G)
        %%
    case 'norne'
        %%
        grdecl = readGRDECL('C:\Users\hnil\Documents\DATA\norne2015\norne\INCLUDE\GRID\IRAP_1005.GRDECL');
        grdecl=cutGrdecl(grdecl,[10 25;35 55;1 22])
        if(opt.vertical)
            grdecl_org=verticalGrdecl(grdecl)
        else
            grdecl_org=grdecl;
        end
        G_org=processGRDECL(grdecl_org);
        if(opt.triangulate)
          faces=unique(G_org.cells.faces(any(bsxfun(@eq,G_org.cells.faces(:,2),[5,6]),2),1));
          %faces=1:G.faces.num;
           G_org=triangulateFaces(G_org,faces');  
        end
        grdecl=padGrdecl(grdecl,[true,true,true],[60 50;40 40;10 10]*3,'relative',true)
        % grdecl=padGrdecl(grdecl,[false,false,true],[20 50;40 40;10 10]/10)
        if(opt.vertical)
            grdecl=verticalGrdecl(grdecl)
        end
        grdecl=refineGrdeclLayers(grdecl,[1 1],opt.ref);
        grdecl=refineGrdeclLayers(grdecl,[grdecl.cartDims(3) grdecl.cartDims(3)],opt.ref);
        G=processGRDECL(grdecl,'Tolerance',opt.gtol);
        %G=triangulateFaces(G,
        G=computeGeometry(G)
        if(opt.triangulate)
        faces=unique(G.cells.faces(any(bsxfun(@eq,G.cells.faces(:,2),[5,6]),2),1));
        %faces=1:G.faces.num;
        G=triangulateFaces(G,faces');
        end
        %checkGrid(G)
        %G=mrstGridWithFullMappings(G);
        %
        figure(),clf,plotGrid(G)
    case 'norne_box'
        %%
        grdecl = readGRDECL('C:\Users\hnil\Documents\DATA\norne2015\norne\INCLUDE\GRID\IRAP_1005.GRDECL');
        grdecl=cutGrdecl(grdecl,[10 25;35 55;1 22])
        if(opt.vertical)
            grdecl_org=verticalGrdecl(grdecl)
        else
            grdecl_org=grdecl;
        end
        G_org=processGRDECL(grdecl_org);
        
        grdecl=padGrdecl(grdecl,[true,true,true],[60 50;40 40;10 10]*3,'relative',true)
        % grdecl=padGrdecl(grdecl,[false,false,true],[20 50;40 40;10 10]/10)
        grdecl=coarseGrdecl(grdecl,[1 1 grdecl.cartDims(end)]);
        if(opt.vertical)
            grdecl=verticalGrdecl(grdecl)
        end
        grdecl=refineGrdecl(grdecl,opt.ref);
        G=processGRDECL(grdecl,'Tolerance',opt.gtol);
        %G=triangulateFaces(G,
        G=computeGeometry(G)
        figure(),clf,plotGrid(G)    
        %%
    otherwise
        error()
end
end