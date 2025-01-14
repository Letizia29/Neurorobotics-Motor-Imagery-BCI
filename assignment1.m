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

% Loading Laplacian Mask
load('laplacian16.mat');

%% Data loading
% for each subject, loading of each run and application of Laplacian filter

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

            % Laplacian masking (pre-processing)
            s = s(:,1:16)*lap;

            % save data
            % offline run
            if run_name(21:27) == 'offline'
                count_off = count_off + 1;
                field_name = strcat("offline",string(count_off));
                subjects.(subj_name).(field_name).s = s;
                subjects.(subj_name).(field_name).h = h;
            end
            % online run
            if run_name(21:26) == 'online'
                count_on = count_on + 1;
                field_name = strcat("online",string(count_on));
                subjects.(subj_name).(field_name).s = s;
                subjects.(subj_name).(field_name).h = h;
            end
        end
    end
end


%% Process the data through log band power computation and ERD/ERS

% for each subject, concatenate the offline files, compute log band power
% and ERD/ERS for each trial, compute the average among trials
% Grand average (whole population analysis): average the values found among
% subjects
% analize if some subjects are far from the avg

sample_rate = 512; % Hz

%% Concatenate signals and events' types, positions, durations for each subject
for i = 1:length(data)
    subj_name = data(i); % one subject
    runs_names = fieldnames(subjects.(subj_name)); % his runs
    subjects.(subj_name).s_c = [];
    subjects.(subj_name).POS = [];
    subjects.(subj_name).DUR = [];
    subjects.(subj_name).TYP = [];

    for j = 1:length(runs_names)
        
        if runs_names{j}(1:7) == 'offline'  
            subjects.(subj_name).DUR = [subjects.(subj_name).DUR; subjects.(subj_name).(runs_names{j}).h.EVENT.DUR];
            subjects.(subj_name).TYP = [subjects.(subj_name).TYP; subjects.(subj_name).(runs_names{j}).h.EVENT.TYP];
            subjects.(subj_name).POS = [subjects.(subj_name).POS; subjects.(subj_name).(runs_names{j}).h.EVENT.POS + length(subjects.(subj_name).s_c)];
            subjects.(subj_name).s_c = [subjects.(subj_name).s_c; subjects.(subj_name).(runs_names{j}).s];
        
        end

    end
    
end

%% Data processing






%% Feature selection
% Identify and extract the most relevant features for each subject and on
% subject average


% % Concatenate the files
% 
% PSD = [file1.PSD; file2.PSD; file3.PSD];
% 
% % Concatenate the events
% h.EVENT.DUR = [file1.h.EVENT.DUR; file2.h.EVENT.DUR; file3.h.EVENT.DUR];
% h.EVENT.TYP = [file1.h.EVENT.TYP; file2.h.EVENT.TYP; file3.h.EVENT.TYP]; 
% h.EVENT.POS = [file1.h.EVENT.POS; file2.h.EVENT.POS + size(file1.PSD, 1); file3.h.EVENT.POS + size(file1.PSD, 1) + size(file2.PSD, 1)];





























