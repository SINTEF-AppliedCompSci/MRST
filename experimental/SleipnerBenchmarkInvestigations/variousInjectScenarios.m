% run several different injection scenarios, each one being written to .mat
% file


ref = [2, 3, 4, 6];

for i=1:numel(ref)

% IEAGHG model grid
studySleipnerBenchmarkFUN('mycase','useIEAGHG_model', 'refineLevel',ref(i));

% studySleipnerBenchmarkFUN('mycase','useIEAGHG_model', 'refineLevel',ref(i), ...
%     'mod_rock_perm', true);
% 
% studySleipnerBenchmarkFUN('mycase','useIEAGHG_model', 'refineLevel',ref(i), ...
%     'mod_rock_poro', true);
% 
% studySleipnerBenchmarkFUN('mycase','useIEAGHG_model', 'refineLevel',ref(i), ...
%     'mod_rhoCO2', true);

studySleipnerBenchmarkFUN('mycase','useIEAGHG_model', 'refineLevel',ref(i), ...
    'mod_rock_perm', true, 'mod_rock_poro', true, 'mod_rhoCO2', true);


% Original model grid
studySleipnerBenchmarkFUN('mycase','useOriginal_model', 'refineLevel',ref(i));

% studySleipnerBenchmarkFUN('mycase','useOriginal_model', 'refineLevel',ref(i), ...
%     'mod_rock_perm', true);
% 
% studySleipnerBenchmarkFUN('mycase','useOriginal_model', 'refineLevel',ref(i), ...
%     'mod_rock_poro', true);
% 
% studySleipnerBenchmarkFUN('mycase','useOriginal_model', 'refineLevel',ref(i), ...
%     'mod_rhoCO2', true);

studySleipnerBenchmarkFUN('mycase','useOriginal_model', 'refineLevel',ref(i), ...
    'mod_rock_perm', true, 'mod_rock_poro', true, 'mod_rhoCO2', true);

end
