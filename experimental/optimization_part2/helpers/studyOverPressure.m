%% Study the over-burden pressure of all formations

% To ensure we do not fracture the caprock during co2 injection, we want to
% keep the bhp of the wells less than the over-burden pressure. Thus, we
% limit the bhp of the wells to be 90 percent of the over-burden pressure,
% as done in Nordbotten & Celia, 2012.

% The over-burden pressure varies with caprock depth, thus we check the max
% and min values within each formation. This range can help identify which
% formations have highly variable caprock depth, and which formations will
% be most limited in terms of storage capacity due to a strict (i.e., low)
% pressure-limit.


% Load formation names:
names = [getBarentsSeaNames() getNorwegianSeaNames() getNorthSeaNames()];

% Remove certain formation names:
names = names(~strcmpi(names,'Nordmelafm'));
names = names(~strcmpi(names,'Rorfm'));
names = names(~strcmpi(names,'Notfm'));
names = names(~strcmpi(names,'Knurrfm'));       % @@ can be constructed??
names = names(~strcmpi(names,'Fruholmenfm'));   % @@
names = names(~strcmpi(names,'Cookfm'));
names = names(~strcmpi(names,'Dunlingp'));
names = names(~strcmpi(names,'Paleocene'));


for i=1:numel(names)
    
    fmName      = names{i};
    fprintf('Processing %s ... \n', fmName);
    coarsening  = 5; % @@ prevent break that occurs when resolution too low.
    rhoCref     = 760 * kilogram / meter ^3;

    
    %%% Get formation
    clear Gt rock2D
    [Gt, rock2D]    = getFormationTopGrid(fmName, coarsening);
    if any(isnan(rock2D.perm))
        rock2D.perm = 500*milli*darcy * ones(Gt.cells.num,1);
    end
    if any(isnan(rock2D.poro))
        rock2D.poro = 0.25 * ones(Gt.cells.num,1); 
    end
    clear seainfo
    seainfo         = getSeaInfo(fmName, rhoCref);

    %%% Get over-burden pressure:
    P_over = computeOverburdenPressure(Gt, rock2D, ...
        seainfo.seafloor_depth, seainfo.water_density);


    % Get max and min over-burden pressure:
    P_over_max(i) = convertTo(max(P_over), mega * Pascal);
    P_over_min(i) = convertTo(min(P_over), mega * Pascal);


    % Compute pressure-limit:
    P_limit_max = P_over_max * 0.9;
    P_limit_min = P_over_min * 0.9;



    % Show in figures:



end

% Show in a table:
fprintf('\n\nOverpressure in MPa:\n')
fprintf('------------------------------------------------\n');
fprintf('\n%-13s&  P-over (max)  &  P-over (min)  \\\\ \n', 'Name');
fprintf('-------------|----------------|----------------\\\\ \n');
for i = 1:numel(names)
    fprintf('%-13s&           %4.0f &           %4.0f \\\\ \n', names{i}, P_over_max(i), P_over_min(i) );
end
fprintf('-------------|----------------|----------------\n\n');

% Show in a bar plot:
figure; set(gcf,'Position',[1 1 1245 400])
bar(1:numel(names), [P_over_min; P_over_max]'); legend('min', 'max')
ylabel('Over-pressure [MPa]', 'FontSize',16)
set(gca,'xtick',1:numel(names),'xticklabel',names,'xticklabelrotation',45, 'fontsize',16)





