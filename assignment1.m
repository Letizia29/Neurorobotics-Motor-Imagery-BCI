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
    POS = [];
    DUR = [];
    TYP = [];

    for j = 1:length(runs_names)
        
        if runs_names{j}(1:7) == 'offline'  
            DUR = [DUR; subjects.(subj_name).(runs_names{j}).h.EVENT.DUR];
            TYP = [TYP; subjects.(subj_name).(runs_names{j}).h.EVENT.TYP];
            POS = [POS; subjects.(subj_name).(runs_names{j}).h.EVENT.POS + length(subjects.(subj_name).s_c)];
            subjects.(subj_name).s_c = [subjects.(subj_name).s_c; subjects.(subj_name).(runs_names{j}).s];
        
        end

    end

    subjects.(subj_name).h.POS = POS;
    subjects.(subj_name).h.DUR = DUR;
    subjects.(subj_name).h.TYP = TYP;

    subjects.(subj_name).vectors = labelVecs(subjects.(subj_name).s_c, subjects.(subj_name).h);
    
end

%% Data processing
% Butterworth filters
n_mu = 4; % filter order
n_beta = 4;

% mu band
W1 = 8; % Hz
W2 = 12; % Hz
Wn_mu = 2*[W1 W2]/sample_rate;

% beta band
W1 = 18; % Hz
W2 = 22; % Hz
Wn_beta = 2*[W1 W2]/sample_rate;

% filter coefficients
[b_mu, a_mu] = butter(n_mu, Wn_mu);
[b_beta, a_beta] = butter(n_beta, Wn_beta);

for i = 1:length(data)
    subj_name = data(i); % one subject
    
    signal = subjects.(subj_name).s_c;
    sfilt_mu = zeros(size(signal));
    sfilt_beta = zeros(size(signal));
    sfilt_sq_mu = zeros(size(signal));
    sfilt_sq_beta = zeros(size(signal));


    for channel = 1:16
        % filter the signal
        sfilt_mu(:, channel) = filtfilt(b_mu, a_mu, signal(:, channel));
        sfilt_beta(:, channel) = filtfilt(b_beta, a_beta, signal(:, channel));
    
        % Rectifying the signal (squaring)
        sfilt_sq_mu(:, channel) = sfilt_mu(:, channel).^2;
        sfilt_sq_beta(:, channel) = sfilt_beta(:, channel).^2; 
    end
    
    % Applying moving average (1-second window)
    LengthWin = 1; % second
    % as FIR filter
    A = 1;
    B = ones(1, LengthWin*h.SampleRate)/LengthWin/h.SampleRate;
    % filter the signal
    sfilt_sq_ma_mu = filter(B, A, sfilt_sq_mu);
    sfilt_sq_ma_beta = filter(B, A, sfilt_sq_beta);
    
    % Logarithm transform
    logBP_mu = log(sfilt_sq_ma_mu);
    logBP_beta = log(sfilt_sq_ma_beta);

    % save results
    subjects.(subj_name).logBP_mu = logBP_mu;
    subjects.(subj_name).logBP_beta = logBP_beta;

end 



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





























