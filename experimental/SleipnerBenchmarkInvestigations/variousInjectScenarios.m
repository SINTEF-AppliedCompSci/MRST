% run several different injection scenarios, each one being written to .mat
% file



studySleipnerBenchmarkFUN('mycase','useOriginal_model', 'refineLevel',-6);

studySleipnerBenchmarkFUN('mycase','useOriginal_model', 'refineLevel',-6, ...
    'mod_rock_perm', true);

studySleipnerBenchmarkFUN('mycase','useOriginal_model', 'refineLevel',-6, ...
    'mod_rock_poro', true);

studySleipnerBenchmarkFUN('mycase','useOriginal_model', 'refineLevel',-6, ...
    'mod_rhoCO2', true);

studySleipnerBenchmarkFUN('mycase','useOriginal_model', 'refineLevel',-6, ...
    'mod_rock_perm', true, 'mod_rock_poro', true, 'mod_rhoCO2', true);