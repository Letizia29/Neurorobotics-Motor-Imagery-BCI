clear
close all
clc

%% Neurorobotics Homework - 1

% Paths

% Custom functions folder path
addpath(genpath('C:\Users\User\Documents\Università\Laurea Magistrale\II Anno\Primo Semestre\Neurorobotics\Homeworks\HW1\functions'))

% biosig package path
addpath(genpath('C:\Users\User\Documents\Università\Laurea Magistrale\Toolboxes\biosig\biosig\t200_FileAccess'))
addpath(genpath('C:\Users\User\Documents\Università\Laurea Magistrale\Toolboxes\biosig\biosig\t250_ArtifactPreProcessingQualityControl'))


%% Load data into a single struct

folder_list = dir('a*'); % Gets the list of all subjects folders

subj = struct();

for i = 1:length(folder_list)
    folder_name = folder_list(i).name;
    subj(i).subj_i = loadSubjData(folder_name);
end


%% Data loading

addpath(genpath('C:\Users\User\Documents\Università\Laurea Magistrale\II Anno\Primo Semestre\Neurorobotics\Homeworks\HW1\HW_data\micontinuous\aj1_micontinuous'));

%% Grand average (whole population analysis)

% Possibile idea:
% Considerare (come ha detto lui) ERD e ERS. Si potrebbero creare degli
% indici per verificare se un soggetto che si discosta di più dalla media
% ha delle feature particolari.

%% DATA LOADING: non viene caricato il file sopen.mat, non so cosa sia

%% DATA SUBJ 1

offline_file1 = 'aj1.20180312.113542.offline.mi.mi_bhbf';
offline_file2 = 'aj1.20180312.114418.offline.mi.mi_bhbf';

online_file1 = 'aj1.20180312.121602.online.mi.mi_bhbf.ema';
online_file2 = 'aj1.20180312.122731.online.mi.mi_bhbf.dynamic';
online_file3 = 'aj1.20180507.112603.online.mi.mi_bhbf.ema';
online_file4 = 'aj1.20180507.113054.online.mi.mi_bhbf.dynamic';
online_file5 = 'aj1.20180601.115626.online.mi.mi_bhbf.ema';
online_file6 = 'aj1.20180601.120114.online.mi.mi_bhbf.dynamic';

% subj1.offline.offline_file1 = offline_file1;
% subj1.offline.offline_file2 = offline_file2;
% 
% subj1.online.online_file1 = online_file1;
% subj1.online.online_file2 = online_file2;
% subj1.online.online_file3 = online_file3;
% subj1.online.online_file4 = online_file4;
% subj1.online.online_file5 = online_file5;
% subj1.online.online_file6 = online_file6;

file_to_load = strcat(offline_file1, '.gdf');

[s, h] = sload(file_to_load);


%% DATA SUBJ 3

offline_file1 = 'aj3.20180313.113110.offline.mi.mi_bhbf';
offline_file2 = 'aj3.20180313.114118.offline.mi.mi_bhbf';
offline_file3 = 'aj3.20180313.114946.offline.mi.mi_bhbf';

online_file1 = 'aj3.20180313.120839.online.mi.mi_bhbf.ema';
online_file2 = 'aj3.20180313.121507.online.mi.mi_bhbf.dynamic';
online_file3 = 'aj3.20180425.111824.online.mi.mi_bhbf.ema';
online_file4 = 'aj3.20180425.112235.online.mi.mi_bhbf.dynamic';
online_file5 = 'aj3.20180529.101640.online.mi.mi_bhbf.ema';
online_file6 = 'aj3.20180529.102142.online.mi.mi_bhbf.dynamic';

%% DATA SUBJ 4

offline_file1 = 'aj4.20180313.151634.offline.mi.mi_bhbf';
offline_file2 = 'aj4.20180313.152525.offline.mi.mi_bhbf';
offline_file3 = 'aj4.20180313.153339.offline.mi.mi_bhbf';

online_file1 = 'aj4.20180313.154650.online.mi.mi_bhbf.ema';
online_file2 = 'aj4.20180313.155139.online.mi.mi_bhbf.dynamic';
online_file3 = 'aj4.20180501.110901.online.mi.mi_bhbf.ema';
online_file4 = 'aj4.20180501.111938.online.mi.mi_bhbf.ema';     % .ema ???? non dovrebbe essere dynamic? chissà
online_file5 = 'aj4.20180522.112515.online.mi.mi_bhbf.ema';
online_file6 = 'aj4.20180522.113943.online.mi.mi_bhbf.dynamic';


%% DATA SUBJ 6

offline_file1 = 'ai6.20180316.153104.offline.mi.mi_bhbf';
offline_file2 = 'ai6.20180316.154006.offline.mi.mi_bhbf';
offline_file3 = 'ai6.20180316.154811.offline.mi.mi_bhbf';

online_file1 = 'ai6.20180316.160351.online.mi.mi_bhbf.ema';
online_file2 = 'ai6.20180316.161026.online.mi.mi_bhbf.dynamic';
online_file3 = 'ai6.20180417.164812.online.mi.mi_bhbf.ema';
online_file4 = 'ai6.20180417.165259.online.mi.mi_bhbf.dynamic';
online_file5 = 'ai6.20180529.151753.online.mi.mi_bhbf.ema';
online_file6 = 'ai6.20180529.152240.online.mi.mi_bhbf.dynamic';


%% DATA SUBJ 7 - prima sessione

offline_file1 = 'ai7.20180316.102257.offline.mi.mi_bhbf';
offline_file2 = 'ai7.20180316.103209.offline.mi.mi_bhbf';
offline_file3 = 'ai7.20180316.104023.offline.mi.mi_bhbf';

online_file1 = 'ai7.20180316.110108.online.mi.mi_bhbf.ema';
online_file2 = 'ai7.20180316.110748.online.mi.mi_bhbf.dynamic';
online_file3 = 'ai7.20180420.101528.online.mi.mi_bhbf.ema';
online_file4 = 'ai7.20180420.101934.online.mi.mi_bhbf.dynamic';
online_file5 = 'ai7.20180518.101353.online.mi.mi_bhbf.ema';
online_file6 = 'ai7.20180518.101711.online.mi.mi_bhbf.dynamic';


%% DATA SUBJ 7 - seconda sessione

offline_file1 = 'aj7.20180323.161608.offline.mi.mi_bhbf';
offline_file2 = 'aj7.20180323.162629.offline.mi.mi_bhbf';

online_file1 = 'aj7.20180323.163958.online.mi.mi_bhbf.ema';
online_file2 = 'aj7.20180323.165113.online.mi.mi_bhbf.ema';
online_file3 = 'aj7.20180420.170712.online.mi.mi_bhbf.ema';
online_file4 = 'aj7.20180420.171150.online.mi.mi_bhbf.dynamic';
online_file5 = 'aj7.20180518.154751.online.mi.mi_bhbf.ema';
online_file6 = 'aj7.20180518.155307.online.mi.mi_bhbf.dynamic';


%% DATA SUBJ 8

offline_file1 = 'ai8.20180320.112744.offline.mi.mi_bhbf';
offline_file2 = 'ai8.20180320.113734.offline.mi.mi_bhbf';
offline_file3 = 'ai8.20180320.114543.offline.mi.mi_bhbf';

online_file1 = 'ai8.20180320.120206.online.mi.mi_bhbf.ema';
online_file2 = 'ai8.20180320.122511.online.mi.mi_bhbf.dynamic';
online_file3 = 'ai8.20180427.152806.online.mi.mi_bhbf.ema';
online_file4 = 'ai8.20180427.154229.online.mi.mi_bhbf.dynamic';
online_file5 = 'ai8.20180601.153831.online.mi.mi_bhbf.ema';
online_file6 = 'ai8.20180601.154427.online.mi.mi_bhbf.dynamic';


%% DATA SUBJ 9

offline_file1 = 'aj9.20180326.153615.offline.mi.mi_bhbf';
offline_file2 = 'aj9.20180326.154532.offline.mi.mi_bhbf';
offline_file3 = 'aj9.20180326.155323.offline.mi.mi_bhbf';

online_file1 = 'aj9.20180326.160922.online.mi.mi_bhbf.ema';
online_file2 = 'aj9.20180326.161314.online.mi.mi_bhbf.dynamic';
online_file3 = 'aj9.20180430.152635.online.mi.mi_bhbf.ema';
online_file4 = 'aj9.20180430.153028.online.mi.mi_bhbf.dynamic';
online_file5 = 'aj9.20180528.112241.online.mi.mi_bhbf.ema';
online_file6 = 'aj9.20180528.112738.online.mi.mi_bhbf.dynamic';






%% Analysis on the whole population

%% Analysis on subject 1


offline_file1 = 'aj1.20180312.113542.offline.mi.mi_bhbf';
offline_file2 = 'aj1.20180312.114418.offline.mi.mi_bhbf';

% online_file1 = 'aj1.20180312.121602.online.mi.mi_bhbf.ema';
% online_file2 = 'aj1.20180312.122731.online.mi.mi_bhbf.dynamic';
% online_file3 = 'aj1.20180507.112603.online.mi.mi_bhbf.ema';
% online_file4 = 'aj1.20180507.113054.online.mi.mi_bhbf.dynamic';
% online_file5 = 'aj1.20180601.115626.online.mi.mi_bhbf.ema';
% online_file6 = 'aj1.20180601.120114.online.mi.mi_bhbf.dynamic';

% Load offline files

file_to_load = strcat(offline_file1, '.gdf');
[s1, h1] = sload(file_to_load);

file_to_load = strcat(offline_file2, '.gdf');
[s2, h2] = sload(file_to_load);

final_filename = strcat(offline_file1, '.mat');


%% Preprocessing

% Laplacian masking
mask = load('HW_data/laplacian16.mat');

% Remove additional channel
if size(s1, 2) > 16
    s1(:, 17:end) = [];
end

if size(s2, 2) > 16
    s2(:, 17:end) = [];
end


%% Initialize useful variables

fs = h1.SampleRate;

N1 = size(s1, 1);
N2 = size(s2, 1);

nchs = size(s1, 2);

%% Concatenate the files

PSD = [file1.PSD; file2.PSD; file3.PSD];

% Concatenate the events
h.EVENT.DUR = [file1.h.EVENT.DUR; file2.h.EVENT.DUR; file3.h.EVENT.DUR];
h.EVENT.TYP = [file1.h.EVENT.TYP; file2.h.EVENT.TYP; file3.h.EVENT.TYP]; 
h.EVENT.POS = [file1.h.EVENT.POS; file2.h.EVENT.POS + size(file1.PSD, 1); file3.h.EVENT.POS + size(file1.PSD, 1) + size(file2.PSD, 1)];