%% Neurorobotics - Assignment 1 - 2024/2025
% Group 8: Neuronauts
% De Rivo Valentino
% Rossato Letizia

clear all
close all
clc


% Paths

% Custom functions folder path
addpath(fullfile(pwd, "functions/"))

% biosig package path
addpath(fullfile(pwd, "Toolboxes\biosig\biosig\t200_FileAccess/"))
addpath(fullfile(pwd, "Toolboxes\biosig\biosig\t250_ArtifactPreProcessingQualityControl/"))


%% Data loading
% for each subject, we have to load each run and compute their PSD
% maybe we should write a function for PSD computation and saving

addpath(genpath(fullfile(pwd, "Data/")))
m = "_micontinuous";
data = [strcat("aj1",m),strcat("aj3",m),strcat("aj4",m),strcat("aj7",m),strcat("aj9",m),strcat("ai6",m),strcat("ai7",m),strcat("ai8",m)];

subjects = struct(); % where s and h data will be saved for each subject

for i = 1:length(data)
    subj_name = data(i);

    % load runs
    runs = dir(fullfile(pwd, strcat("Data/", subj_name)));
    runs_names = {runs.name};
    
    count_off = 0;
    count_on = 0;
    for j = 1:length(runs_names)
        run_name = runs_names{j};
        if run_name(1) == 'a' % actual run
            % load data
            [s, h] = sload(run_name);
            % save data
            if run_name(21:27) == 'offline'
                count_off = count_off + 1;
                field_name = strcat("offline",string(count_off));
                subjects.(subj_name).(field_name).s = s;
                subjects.(subj_name).(field_name).h = h;
            end

            if run_name(21:26) == 'online'
                count_on = count_on + 1;
                field_name = strcat("online",string(count_on));
                subjects.(subj_name).(field_name).s = s;
                subjects.(subj_name).(field_name).h = h;
            end
        end
    end
end

%% DA AGGIUNGERE DIRETTAMENTE DENTRO IL CICLO?
%% Preprocessing

% Laplacian masking
mask = load('laplacian16.mat');

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


%% Grand average (whole population analysis)

% Possibile idea:
% Considerare (come ha detto lui) ERD e ERS. Si potrebbero creare degli
% indici per verificare se un soggetto che si discosta di più dalla media
% ha delle feature particolari.






























