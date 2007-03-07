%% pack_se
% script to pack up sePARSE code.
% se10
% |- plotter
% |- make_nii, from mathworks user site.
% |- save_nii, from mathworks user site.
% |- dataset_load
% |- se_opts
% |- varparser
% |- load_procpar
% '- query_procpar
% 
clear
d='sePARSE';
mkdir(d)

F{1}=which('plotter');
F{2}=which('dataset_load');
F{3}=which('se_opts');
F{4}=which('varparser');
F{5}=which('load_procpar');
F{6}=which('query_procpar');
F{7}=which('se10');
F{8}=which('pack_se');

for f=F
    s=['cp ' cell2mat(f) ' ' d '/'];
    system(s);
end
s=['tar cvzf sePARSE_10_' date '.tgz sePARSE'];
system(s);
s=['scp sePARSE_10_' date '.tgz markbold@foo.vsrc.uab.edu:~/Sites/sePARSE/'];
system(s);

mlint -cyc se10 % code complexity metric
