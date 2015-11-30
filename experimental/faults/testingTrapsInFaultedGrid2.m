% Calling makeInternalBoundaries() takes the specified faces and makes two
% new faces that are not connected. This causes the indexing of the faces
% to be modified. Rather than using this function, simply modify the
% transmissibilities of the fault faces (0 for sealing, between 0 and 1 for
% conducting). The cells adjacent to the fault faces should ideally have an
% updated elevation, but this is probably only important if faults were
% considered to be leaking. The atlas model does not contain the right info
% to update cell elevation in the same manner as done here? (need to confirm
% this)

moduleCheck('ad-core')

%% Load three different models of Johansen
% (taken from trapsJohansen.m)

% can any faults be loaded into the grdecl?

%%% Load full-field model
[ grdecl_full, G_full, Gt_full, rock2D_full ] = makeJohansenFullField();

%%% Load sector model, and construct rock2D from full-field model
[ grdecl_sect, G_sect, Gt_sect, rock2D_sect ] = makeJohansenSector();
% note that if topSurfaceGrid(...'AddFaults',true) is called, Gt.faces.tag
% appear. However, tags appear for internal faces along outside of domain.
% This should be de-bugged and tested more.

% Get the position of the well (data given from 'Sector5_Well.txt');
wi    = find(G_sect.cells.indexMap==sub2ind(G_sect.cartDims, 48, 48, 1));
wcent = G_sect.cells.centroids(wi,:);
d = sqrt(sum(bsxfun(@minus, Gt_sect.cells.centroids, wcent(1:2)).^2, 2));
[~,wi_sect] = min(d);

%%% Load atlas model
[ grdecl_atls, G_atls, Gt_atls, rock2D_atls ] = makeJohansenAtlas();
[ ~, fault_finx_atls ] = implementFaults(Gt_atls, rock2D_atls,'gradTol',100);
Gt_atls = makeInternalBoundary(Gt_atls, fault_finx_atls); % @@
Gt_atls = computeGeometryVE(Gt_atls);

%% Compare which 3D / 2D grids contain fault info

figure; set(gcf,'name','Full-field','Position',[2562 860 1500 500])
subplot(1,2,1); title('3D grid')
plotGrid(G_full,'EdgeAlpha',0.1,'FaceColor','none'); view(3)
subplot(1,2,2); title('2D grid')
plotGrid(Gt_full, 'EdgeAlpha',0.1,'FaceColor','none'); view(3)

figure; set(gcf,'name','Sector','Position',[3076 438 1500 500])
subplot(1,2,1); title('3D grid')
plotGrid(G_sect,'EdgeAlpha',0.1,'FaceColor','none'); view(3)
subplot(1,2,2); title('2D grid')
plotGrid(Gt_sect, 'EdgeAlpha',0.1,'FaceColor','none'); view(3)

figure; set(gcf,'name','Atlas','Position',[3603 56 1500 500])
subplot(1,2,1); title('3D grid')
plotGrid(G_atls,'EdgeAlpha',0.1,'FaceColor','none'); view(3)
subplot(1,2,2); title('2D grid')
plotGrid(Gt_atls, 'EdgeAlpha',0.1,'FaceColor','none'); view(3)


%% Previous work done on Sector model to add fault to top grid
% (from makeJohansenVEgridSIAM.m)
% Fault cells are manually specified, used to get fault faces, and then
% used to introduce internal boundaries. After internal bdrys added, cell
% elevations along fault must be updated.

% initial cell elevations
figure
plotGrid(Gt_sect, 'FaceColor','none')
plotCellData(Gt_sect, Gt_sect.cells.z, 'EdgeAlpha',0.1)
axis tight off; view(3)
title('Before internal boundaries created')

%%% FIND FAULT - set inner boundary
% The main fault is assumed to be sealing and must therefore be represented
% as an inner boundary in the 2D grid.  Start by locating cells on each
% side of fault, which is defined to be between cells with index i=43 and
% i=44 for j<44.
cells2D_1 = find(Gt_sect.cells.ij(:,1) == 43 & Gt_sect.cells.ij(:,2) <= 44);
cells2D_2 = find(Gt_sect.cells.ij(:,1) == 44 & Gt_sect.cells.ij(:,2) <= 44);

% Plot the cells on opposite sides
figure;
plotGrid(Gt_sect, 'faceColor', 'none');
plotGrid(Gt_sect, cells2D_1, 'faceColor', 'r')
plotGrid(Gt_sect, cells2D_2, 'faceColor', 'g')
axis tight off; view(3)
title('Cells on opposite sides of sealing fault'),

% Find the faces at the fault. Construct a mapping facesMat defined such
% that facesMat(i,j)=k if face <k> is shared by cells <i> and <j>.
facesMat = sparse(double(Gt_sect.faces.neighbors(:,1))+1, ...
   double(Gt_sect.faces.neighbors(:,2))+1, 1:Gt_sect.faces.num);
facesMat = facesMat + facesMat';
faultFaces2D = diag( facesMat(cells2D_1+1, cells2D_2+1) );

% Make internal boundary and compute the geometry of the resulting grid
%Gt_sect = makeInternalBoundary(Gt_sect, faultFaces2D);
%Gt_sect = computeGeometryVE(Gt_sect);

% NB: makeInternalBoundary sets Gt.faces.tag = -1 for any internal bdry.
% Get fault_finx using tags, and plot fault faces.
%fault_finx = [1:Gt_sect.faces.num]';
%fault_finx = fault_finx( Gt_sect.faces.tag<0 );

figure
plotGrid(Gt_sect, 'faceColor', 'none');
plotCellData(Gt_sect, Gt_sect.cells.z, 'EdgeAlpha',0.1)
plotFaces(Gt_sect, fault_finx, 'EdgeColor','r','LineWidth',3) % @@
axis tight off; view(3)
title('After internal boundaries created: incorrect Gt.cells.z')

% The function 'topSurfaceGrid' does not handle faults correctly when the
% cells on opposite sides are not phyiscally in contact with each other.
% Instead of producing a discontinuity in the 'z' value, an average value
% is used.  Hence, we need to manually reset the 'z' value of these cells
% (marked in red and green in the plot) to avoid an incorrect flow 
Gt_sect.cells.z([cells2D_1; cells2D_2]) = ...
   G_sect.faces.centroids(Gt_sect.cells.map3DFace([cells2D_1; cells2D_2]), 3);

figure
plotGrid(Gt_sect, 'faceColor', 'none');
plotCellData(Gt_sect, Gt_sect.cells.z, 'EdgeAlpha',0.1)
plotFaces(Gt_sect, fault_finx, 'EdgeColor','r','LineWidth',3)
axis tight off; view(3)
title('After internal boundaries created, with updated Gt.cells.z')

% NB: other grid data should be modified (like nodes?) in order to
% correctly plot the faulted cells as non-connecting

%% Adjust transmissibilities
s = setupOperatorsTPFA(Gt_sect, rock2D_sect);

% 
Gt = Gt_sect;
rock = rock2D_sect;

rock_tmp      = rock; 
rock_tmp.perm = rock.perm .* Gt.cells.H; 
T             = computeTrans(Gt, rock_tmp); 
cf            = Gt.cells.faces(:, 1); 
nf            = Gt.faces.num; 
T             = 1 ./ accumarray(cf, 1 ./ T, [nf, 1]);

s.T_all = T;
transMult = 0; %less than 1 if conducting (must be between 0 and 1)
s_fault = s;
%s_fault.T_all( fault_finx ) = s_fault.T_all( fault_finx ).*transMult;
s_fault.T_all( faultFaces2D ) = s_fault.T_all(faultFaces2D ).*transMult;
s_fault.T = s_fault.T_all( s.internalConn );


% should s.internalConn be updated?
% should node elevation be updated? neighbors?

%% Simulate

% Using grid with internal bdrys and adjusted cell elevations
[ initState, fluid ] = simpleSetUp( Gt_sect, rock2D_sect );

% Using adjusted transmissibilities (T=T*transMult for all fault faces)
model = CO2VEBlackOilTypeModel_trans(Gt_sect, rock2D_sect, fluid, 'trans', s_fault.T_all); 

% Using un-modified transmissibilities
model = CO2VEBlackOilTypeModel(Gt_sect, rock2D_sect, fluid);

% Using un-modified grid (fault represented as dramatic incline)
[ grdecl_sect, G_sect, Gt_sect, rock2D_sect ] = makeJohansenSector();
[ initState, fluid ] = simpleSetUp( Gt_sect, rock2D_sect );
model = CO2VEBlackOilTypeModel(Gt_sect, rock2D_sect, fluid);
                                     

itime   = 10*year;
istep   = 10;
qtot    = 1*1e9*kilogram*convertTo(itime,year)/(rhoc*kilogram/meter^3); % total m3 over entire injection time
mtime   = 2*year;
mstep   = 2;

schedule = setSchedule(Gt_sect, rock2D_sect, wi_sect, ...
    qtot, istep, itime, mstep, mtime, true);    % @@ bug if mtime = 1 but mstep = 2;

% Set up boundary conditions
bdryFaces   = find( Gt_sect.faces.neighbors(:,1).*Gt_sect.faces.neighbors(:,2) == 0 );
bdryVal     = Gt_sect.faces.z(bdryFaces) * rhow*kilogram/meter^3 * norm(gravity);
bc          = addBC( [], bdryFaces, 'pressure', bdryVal, 'sat', [1 0] );
for i = 1:numel(schedule.control)
    schedule.control(i).bc = bc;
end


[wellSols, states, sim_report] = simulateScheduleAD(initState, model, schedule);
    

%% Other method to implement fault

[ ~, faultFaces ] = implementFaults(Gt_f, rock2D, 'gradTol',100);

figure
tmp_z = zeros(Gtf.faces.num,1);
tmp_z( faultFaces ) = Gtf.faces.z( faultFaces ); 
plotFaceData(Gtf, tmp_z) % didn't work

figure
plotFaces(Gtf, faultFaces, 'r')
plotGrid(Gtf, 'FaceColor','none','EdgeAlpha',0.1)



% tips:
% - focus on modifying transmissibilities and then running a simulation
% - when creating model, don't pass in rock2D, but set s.T
% - when plotting, plot node.z (after first modifying node.z at any fault
% faces
% - fault could be conducting or sealing, even if fault cells on either
% side of fault are not physically touching


%%
[Gt, rock2D] = getFormationTopGrid('Garnfm',5);
[~, faultFaces, faultCells, s ] = implementFaults(Gt, rock2D);

% pass in s to create model (not rock2D).







