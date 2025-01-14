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

                % Calculate PSD and save relevant information
                [PSD, h_PSD, f] = get_PSD(s, h);

                subjects.(subj_name).(field_name).PSD = PSD;
                subjects.(subj_name).(field_name).h_PSD = h_PSD;
                subjects.(subj_name).(field_name).f = f;
            end
            % online run
            if run_name(21:26) == 'online'
                count_on = count_on + 1;
                field_name = strcat("online",string(count_on));
                subjects.(subj_name).(field_name).s = s;
                subjects.(subj_name).(field_name).h = h;

                % Calculate PSD and save relevant information
                [PSD, h_PSD, f] = get_PSD(s, h);

                subjects.(subj_name).(field_name).PSD = PSD;
                subjects.(subj_name).(field_name).h_PSD = h_PSD;
                subjects.(subj_name).(field_name).f = f;
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

%% ERD/ERS Spectrogram

%% Concatenate the files

for i = 1:length(data)
    subj_name = data(i);
    runs_names = fieldnames(subjects.(subj_name));
    subjects.(subj_name).PSD_c = [];
    POS = [];
    DUR = [];
    PSD = [];

    for j = 1:length(runs_names)

        if contains(runs_names{j}, 'offline', 'IgnoreCase', true)
        
        % QUESTO NON FUNZIONA PIU' PERCHE' I FILE CONTENGONO NOMI LUNGHI 3
        % CARATTERI
        %if runs_names{j}(1:7) == 'offline'  

            DUR = [DUR; subjects.(subj_name).(runs_names{j}).h_PSD.EVENT.DUR];
            POS = [POS; subjects.(subj_name).(runs_names{j}).h_PSD.EVENT.POS + length(subjects.(subj_name).PSD_c)];

            subjects.(subj_name).PSD_c = [subjects.(subj_name).PSD_c; subjects.(subj_name).(runs_names{j}).PSD];        
        end

    end

    subjects.(subj_name).h_PSD.POS = POS;
    subjects.(subj_name).h_PSD.DUR = DUR;
    subjects.(subj_name).h_PSD.TYP = subjects.(subj_name).h.TYP;

    subjects.(subj_name).vectors_PSD = labelVecs(subjects.(subj_name).PSD_c, subjects.(subj_name).h_PSD);
    
end

%% Activity and Reference computation

ax1 = figure(1); 
% ax1.Name = '';
ax1.Position = [50,100,1600,600];

for i = 1:length(data)
    subj_name = data(i);
    runs_names = fieldnames(subjects.(subj_name));

    subjects.(subj_name).ERD = [];

    % Store the values in temp variables
    h_PSD = subjects.(subj_name).h_PSD;
    PSD_c = subjects.(subj_name).PSD_c;

    % Get the starting position of the trials
    startTrial = h_PSD.POS(h_PSD.TYP == 786);
    stopTrial  = h_PSD.POS(h_PSD.TYP == 781) + h_PSD.DUR(h_PSD.TYP == 781);

    % Get the number and length of trials
    ntrials = length(startTrial);
    trial_length = min(stopTrial - startTrial);

    Activity = zeros(trial_length, size(PSD_c, 2), size(PSD_c, 3), ntrials);      % [windows x frequencies x channels x trials]

    for trId = 1 : ntrials
        cstart = startTrial(trId);
        cstop = stopTrial(trId);

        % PSD
        Activity(:, :, :, trId) = PSD_c(cstart:min(cstop, cstart+trial_length-1), :, :);
    end

    % Extract the data referring to the fixation period
    durFixData  = min(h_PSD.DUR(h_PSD.TYP == 786));

    FixData = Activity(1:durFixData, :, :, :);

    % Reference
    Reference = repmat(mean(FixData), [size(Activity, 1) 1 1 1]);

    %subjects.(subj_name).ERD = 100 * (Activity- Reference)./ Reference;
    ERD = log(Activity ./ Reference);
    subjects.(subj_name).ERD = ERD;


    Ck = h_PSD.TYP(h_PSD.TYP == 771 | h_PSD.TYP == 773);

    ERDavg_feet  = mean(ERD(:, :, :, Ck == 771), 4);
    ERDavg_hands = mean(ERD(:, :, :, Ck == 773), 4);


    % Visualization

    chns = [7, 9, 11];

    f = subjects.(subj_name).(runs_names{1}).f;

    ERDavg_feet_rot = imrotate(ERDavg_feet(:, :, 9), 90);

    % FARE QUALCOSA SU QUESTE OSCENITA'
    subplot(2, 4, mod(i-1, 8)+1) ;
    imagesc(0:8, f, ERDavg_feet_rot);
    colormap hot
    colorbar
    clim([-1.1, 1.7])
    %xlim([0, 48])
    %ylim([0, 50])
    set(gca,'YDir','reverse')
    
end


%% Feature selection

% Identify and extract the most relevant features for each subject and on
% subject average

