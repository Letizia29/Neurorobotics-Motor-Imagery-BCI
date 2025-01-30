%% Neurorobotics - Assignment 1 - 2024/2025
% Group 8: Neuronauts
% De Rivo Valentino
% Rossato Letizia

clear
close all
clc

% Paths

% Custom functions folder path
addpath(fullfile(pwd, "functions/"))

% biosig package path
addpath(fullfile(pwd, "Toolboxes\biosig\biosig\t200_FileAccess/"))
addpath(fullfile(pwd, "Toolboxes\biosig\biosig\t250_ArtifactPreProcessingQualityControl/"))

% Data path
addpath(genpath(fullfile(pwd, "Data/")))

% Add custom path to EEGLAB


%% Data loading and PSD computation

% Define channels' labels
channels = {"Fz", "FC3", "FC1", "FCz", "FC2", "FC4", "C3", "C1", "Cz", "C2", "C4", "CP3", "CP1", "CPz", "CP2", "CP4"};

% Loading Laplacian Mask
load('laplacian16.mat');

% For each subject, loading of each run, application of Laplacian filter,
% artifacts removal and PSD computation, events re-computation

m = "_micontinuous";
data = [strcat("aj1",m),strcat("aj3",m),strcat("aj4",m),strcat("aj7",m),strcat("aj9",m),strcat("ai6",m),strcat("ai7",m),strcat("ai8",m)];

% Initialize struct to save subjects' data
subjects = struct(); 

disp('Loading subjects data')

% Count the number of removed artifacts
tot_num_rem_trials = 0;

for i = 1:length(data)
    subj_name = data(i);

    % Load runs
    runs = dir(fullfile(pwd, strcat("Data/", subj_name)));
    runs_names = {runs.name};
    
    count_off = 0;
    count_on = 0;
    for j = 1:length(runs_names)
        run_name = runs_names{j};
        if run_name(1) == 'a'       % Run

            % Load data
            [s, h] = sload(run_name);

            % Laplacian masking (pre-processing)
            s = s(:,1:16)*lap;

            % Save data
            % offline run
            if run_name(21:27) == 'offline'
                count_off = count_off + 1;
                field_name = strcat("offline",string(count_off));

                % Signal filtering
                Wn1 = 3/(512/2);
                Wn2 = 49/(512/2);
                [b, a] = butter(3, [Wn1, Wn2]);
                s = filtfilt(b, a, s);

                % Remove artifacts
                [s, h, numremtrials] = remArtifacts(s, h);

                % Update counter
                tot_num_rem_trials = tot_num_rem_trials + numremtrials;

                disp(['N_trials removed from subj ', data{i}(1:3), ' run ', num2str(j-2), ': ', num2str(numremtrials)])

                % Save information for each run
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

                % Save information for each run
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

disp(['Number of trials removed in total: ', num2str(tot_num_rem_trials)])

%% Process the data through log band power computation and ERD/ERS

% for each subject, concatenate the offline files, compute log band power
% and ERD/ERS for each trial, compute the average among trials
% Grand average (whole population analysis): average the values found among
% subjects
% analize if some subjects are far from the avg

sample_rate = 512; % [Hz]


%% Concatenate signals and events' types, positions, durations for each subject

clc

for i = 1:length(data)

    disp(['Subject ', data{i}(1:3)])

    subj_name = data(i);                                % Subject
    runs_names = fieldnames(subjects.(subj_name));      % Runs

    % Initialize structure fields
    subjects.(subj_name).s_c = [];
    POS = [];
    DUR = [];
    TYP = [];

    for j = 1:length(runs_names)
        
        if runs_names{j}(1:7) == 'offline'

            % Concatenation of the offline runs for each subject
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

clc
disp('Done')

%% Data processing

% Butterworth filters
n_mu = 4;               % Filter order
n_beta = 4;

% mu band
W1 = 8;                 % [Hz]
W2 = 12;                % [Hz]
Wn_mu = 2*[W1 W2]/sample_rate;

% beta band
W1 = 16;                % [Hz]
W2 = 24;                % [Hz]
Wn_beta = 2*[W1 W2]/sample_rate;

% Filter coefficients
[b_mu, a_mu] = butter(n_mu, Wn_mu);
[b_beta, a_beta] = butter(n_beta, Wn_beta);

% Store starting time information for Cue and Continuous Feedback
starting_time_cue = zeros(length(data), 1);
starting_time_cf  = zeros(length(data), 1);

for i = 1:length(data)

    disp(['Processing subject ', data{i}(1:3), ' data'])

    subj_name = data(i);
    
    % Initialize data structures
    signal = subjects.(subj_name).s_c;
    sfilt_mu = zeros(size(signal));
    sfilt_beta = zeros(size(signal));
    sfilt_sq_mu = zeros(size(signal));
    sfilt_sq_beta = zeros(size(signal));

    for channel = 1:16
        % Filter the signal
        sfilt_mu(:, channel)   = filtfilt(b_mu, a_mu, signal(:, channel));
        sfilt_beta(:, channel) = filtfilt(b_beta, a_beta, signal(:, channel));
    
        % Rectify the signal (squaring)
        sfilt_sq_mu(:, channel)   = sfilt_mu(:, channel).^2;
        sfilt_sq_beta(:, channel) = sfilt_beta(:, channel).^2; 
    end
    
    % Applying moving average (1-second window)
    LengthWin = 1;          % [s]
    % as FIR filter
    A = 1;
    B = ones(1, LengthWin*sample_rate)/LengthWin/sample_rate;
    % Apply the filter
    sfilt_sq_ma_mu = filter(B, A, sfilt_sq_mu);
    sfilt_sq_ma_beta = filter(B, A, sfilt_sq_beta);
    
    % LOG BAND POWER
    % Logarithm transform
    logBP_mu = log(sfilt_sq_ma_mu);
    logBP_beta = log(sfilt_sq_ma_beta);

    % Save results
    subjects.(subj_name).logBP_mu = logBP_mu;
    subjects.(subj_name).logBP_beta = logBP_beta;
    
    % ERD/ERS ON LOGARITHMIC BAND POWER
    % Extract trials for the 2 classes
    % Get the starting position of the trials
    h = subjects.(subj_name).h;
    startTrial = h.POS(h.TYP == 786);
    stopTrial  = h.POS(h.TYP == 781) + h.DUR(h.TYP == 781);

    % Get the number and length of trials
    ntrials = length(startTrial);
    trial_length = min(stopTrial - startTrial);

    % Find the minimum length trial
    min_trial = find(trial_length == (stopTrial - startTrial), 1, "first");

    % Save the starting time of Cue and Continuous Feedback
    starting_time_cue(i) = h.DUR(min_trial*4-2);
    starting_time_cf(i)  = h.DUR(min_trial*4-2) + h.DUR(min_trial*4-1);

    Activity_mu = zeros(trial_length, size(sfilt_sq_ma_mu, 2), ntrials);          % [samples x channels x trials]
    Activity_beta = zeros(trial_length, size(sfilt_sq_ma_beta, 2), ntrials);      % [samples x channels x trials]

    for trId = 1 : ntrials
        cstart = startTrial(trId);
        cstop = stopTrial(trId);

        % Activity periods
        Activity_mu(:, :, trId) = sfilt_sq_ma_mu(cstart:min(cstop, cstart+trial_length-1), :);
        Activity_beta(:, :, trId) = sfilt_sq_ma_beta(cstart:min(cstop, cstart+trial_length-1), :);
    end

    % Extract the data referring to the fixation period
    durFixData  = min(h.DUR(h.TYP == 786));

    FixData_mu = Activity_mu(1:durFixData, :, :);
    FixData_beta = Activity_beta(1:durFixData, :, :);

    % Reference
    Reference_mu = repmat(mean(FixData_mu), [size(Activity_mu, 1) 1 1]);
    Reference_beta = repmat(mean(FixData_beta), [size(Activity_beta, 1) 1 1]);

    % Calculate and store ERD_logBP values
    subjects.(subj_name).ERD_logBP_mu = 100 * (Activity_mu - Reference_mu) ./ Reference_mu;
    subjects.(subj_name).ERD_logBP_beta = 100 * (Activity_beta - Reference_beta) ./ Reference_beta;

end 

clc
disp('Done')

%% Plot ERD/ERS on logBP

% Significant channels (C3, Cz, C4)
chns = [7, 9, 11];

ax1 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax1.Name = 'Average ERD logBP \mu - channel C3'; % 7
ax1.NumberTitle = 'off';
hold off

ax2 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax2.Name = 'Average ERD logBP \mu - channel Cz'; % 9
ax2.NumberTitle = 'off';
hold off

ax3 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax3.Name = 'Average ERD logBP \mu - channel C4'; % 11
ax3.NumberTitle = 'off';
hold off

ax4 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax4.Name = 'Average ERD logBP \beta - channel C3'; % 7
ax4.NumberTitle = 'off';
hold off

ax5 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax5.Name = 'Average ERD logBP \beta - channel Cz'; % 9
ax5.NumberTitle = 'off';
hold off

ax6 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax6.Name = 'Average ERD logBP \beta - channel C4'; % 11
ax6.NumberTitle = 'off';
hold off

for i = 1:length(data)
    subj_name = data(i);

    % Separate the classes for each frequency band
    ERDmu_feet    = subjects.(subj_name).ERD_logBP_mu(:, :, subjects.(subj_name).vectors.Ck == 771);
    ERDmu_hands   = subjects.(subj_name).ERD_logBP_mu(:, :, subjects.(subj_name).vectors.Ck == 773);
    ERDbeta_feet  = subjects.(subj_name).ERD_logBP_beta(:, :, subjects.(subj_name).vectors.Ck == 771);
    ERDbeta_hands = subjects.(subj_name).ERD_logBP_beta(:, :, subjects.(subj_name).vectors.Ck == 773);

    % Average ERD/ERS
    ERD_logBP_mu_avg_feet  = mean(ERDmu_feet, 3);
    ERD_logBP_mu_avg_hands = mean(ERDmu_hands, 3);

    ERD_logBP_beta_avg_feet  = mean(ERDbeta_feet, 3);
    ERD_logBP_beta_avg_hands = mean(ERDbeta_hands, 3);

    % Standard error
    ERD_logBP_mu_SE_feet  = std(ERDmu_feet, 0, 3)./sqrt(length(ERDmu_feet(1,1,:)));
    ERD_logBP_mu_SE_hands = std(ERDmu_hands, 0, 3)./sqrt(length(ERDmu_hands(1,1,:)));

    ERD_logBP_beta_SE_feet  = std(ERDbeta_feet, 0, 3)./sqrt(length(ERDbeta_feet(1,1,:)));
    ERD_logBP_beta_SE_hands = std(ERDbeta_hands, 0, 3)./sqrt(length(ERDbeta_hands(1,1,:)));


    % Visualization

    % Time vector
    T = 1/sample_rate;                              % [s]
    t = 0:T:(length(ERD_logBP_mu_avg_feet)-1)*T;    % [s]

    set(0,'CurrentFigure',ax1)
    subplot(2, 4, mod(i-1, 8)+1)
    plot(t, ERD_logBP_mu_avg_feet(:, chns(1)), 'g'), hold on
    plot(t, ERD_logBP_mu_avg_hands(:, chns(1)), 'r')
    % SE
    plot(t, ERD_logBP_mu_avg_hands(:, chns(1)) - ERD_logBP_mu_SE_hands(:, chns(1)), ':r')
    plot(t, ERD_logBP_mu_avg_hands(:, chns(1)) + ERD_logBP_mu_SE_hands(:, chns(1)), ':r')
    plot(t, ERD_logBP_mu_avg_feet(:, chns(1)) - ERD_logBP_mu_SE_feet(:, chns(1)), ':g')
    plot(t, ERD_logBP_mu_avg_feet(:, chns(1)) + ERD_logBP_mu_SE_feet(:, chns(1)), ':g')
    % Starting time of Cue and CF
    time_cue = starting_time_cue(i)/512;
    time_cf = starting_time_cf(i)/512;
    xline(time_cue, 'k--')
    xline(time_cf, 'k--')
    name = char(subj_name);
    title(strcat('ERD avg logBP \mu - C3 - ',name(1:3)))
    legend('both feet', 'both hands')
    axis tight
    drawnow
    hold off

    set(0,'CurrentFigure',ax2)
    subplot(2, 4, mod(i-1, 8)+1)
    plot(t, ERD_logBP_mu_avg_feet(:, chns(2)), 'g'), hold on
    plot(t, ERD_logBP_mu_avg_hands(:, chns(2)), 'r')
    % SE
    plot(t, ERD_logBP_mu_avg_hands(:, chns(2)) - ERD_logBP_mu_SE_hands(:, chns(2)), ':r')
    plot(t, ERD_logBP_mu_avg_hands(:, chns(2)) + ERD_logBP_mu_SE_hands(:, chns(2)), ':r')
    plot(t, ERD_logBP_mu_avg_feet(:, chns(2)) - ERD_logBP_mu_SE_feet(:, chns(2)), ':g')
    plot(t, ERD_logBP_mu_avg_feet(:, chns(2)) + ERD_logBP_mu_SE_feet(:, chns(2)), ':g')
    % Starting time of Cue and CF
    time_cue = starting_time_cue(i)/512;
    time_cf = starting_time_cf(i)/512;
    xline(time_cue, 'k--')
    xline(time_cf, 'k--')
    name = char(subj_name);
    title(strcat('ERD avg logBP \mu - Cz - ',name(1:3)))
    legend('both feet', 'both hands')
    axis tight
    drawnow
    hold off

    set(0,'CurrentFigure',ax3)
    subplot(2, 4, mod(i-1, 8)+1)
    plot(t, ERD_logBP_mu_avg_feet(:, chns(3)), 'g'), hold on
    plot(t, ERD_logBP_mu_avg_hands(:, chns(3)), 'r')
    % SE
    plot(t, ERD_logBP_mu_avg_hands(:, chns(3)) - ERD_logBP_mu_SE_hands(:, chns(3)), ':r')
    plot(t, ERD_logBP_mu_avg_hands(:, chns(3)) + ERD_logBP_mu_SE_hands(:, chns(3)), ':r')
    plot(t, ERD_logBP_mu_avg_feet(:, chns(3)) - ERD_logBP_mu_SE_feet(:, chns(3)), ':g')
    plot(t, ERD_logBP_mu_avg_feet(:, chns(3)) + ERD_logBP_mu_SE_feet(:, chns(3)), ':g')
    % Starting time of Cue and CF
    time_cue = starting_time_cue(i)/512;
    time_cf = starting_time_cf(i)/512;
    xline(time_cue, 'k--')
    xline(time_cf, 'k--')
    name = char(subj_name);
    title(strcat('ERD avg logBP \mu - C4 - ',name(1:3)))
    legend('both feet', 'both hands')
    axis tight
    drawnow
    hold off

    set(0,'CurrentFigure',ax4)
    subplot(2, 4, mod(i-1, 8)+1)
    plot(t, ERD_logBP_beta_avg_feet(:, chns(1)), 'g'), hold on
    plot(t, ERD_logBP_beta_avg_hands(:, chns(1)), 'r')
    % SE
    plot(t, ERD_logBP_beta_avg_hands(:, chns(1)) - ERD_logBP_beta_SE_hands(:, chns(1)), ':r')
    plot(t, ERD_logBP_beta_avg_hands(:, chns(1)) + ERD_logBP_beta_SE_hands(:, chns(1)), ':r')
    plot(t, ERD_logBP_beta_avg_feet(:, chns(1)) - ERD_logBP_beta_SE_feet(:, chns(1)), ':g')
    plot(t, ERD_logBP_beta_avg_feet(:, chns(1)) + ERD_logBP_beta_SE_feet(:, chns(1)), ':g')
    % Starting time of Cue and CF
    time_cue = starting_time_cue(i)/512;
    time_cf = starting_time_cf(i)/512;
    xline(time_cue, 'k--')
    xline(time_cf, 'k--')
    name = char(subj_name);
    title(strcat('ERD avg logBP \beta - C3 - ',name(1:3)))
    legend('both feet', 'both hands')
    axis tight
    drawnow
    hold off

    set(0,'CurrentFigure',ax5)
    subplot(2, 4, mod(i-1, 8)+1)
    plot(t, ERD_logBP_beta_avg_feet(:, chns(2)), 'g'), hold on
    plot(t, ERD_logBP_beta_avg_hands(:, chns(2)), 'r')
    % SE
    plot(t, ERD_logBP_beta_avg_hands(:, chns(2)) - ERD_logBP_beta_SE_hands(:, chns(2)), ':r')
    plot(t, ERD_logBP_beta_avg_hands(:, chns(2)) + ERD_logBP_beta_SE_hands(:, chns(2)), ':r')
    plot(t, ERD_logBP_beta_avg_feet(:, chns(2)) - ERD_logBP_beta_SE_feet(:, chns(2)), ':g')
    plot(t, ERD_logBP_beta_avg_feet(:, chns(2)) + ERD_logBP_beta_SE_feet(:, chns(2)), ':g')
    % Starting time of Cue and CF
    time_cue = starting_time_cue(i)/512;
    time_cf = starting_time_cf(i)/512;
    xline(time_cue, 'k--')
    xline(time_cf, 'k--')
    name = char(subj_name);
    title(strcat('ERD avg logBP \beta - Cz - ',name(1:3)))
    legend('both feet', 'both hands')
    axis tight
    drawnow
    hold off

    set(0,'CurrentFigure',ax6)
    subplot(2, 4, mod(i-1, 8)+1)
    plot(t, ERD_logBP_beta_avg_feet(:, chns(3)), 'g'), hold on
    plot(t, ERD_logBP_beta_avg_hands(:, chns(3)), 'r')
    % SE
    plot(t, ERD_logBP_beta_avg_hands(:, chns(3)) - ERD_logBP_beta_SE_hands(:, chns(3)), ':r')
    plot(t, ERD_logBP_beta_avg_hands(:, chns(3)) + ERD_logBP_beta_SE_hands(:, chns(3)), ':r')
    plot(t, ERD_logBP_beta_avg_feet(:, chns(3)) - ERD_logBP_beta_SE_feet(:, chns(3)), ':g')
    plot(t, ERD_logBP_beta_avg_feet(:, chns(3)) + ERD_logBP_beta_SE_feet(:, chns(3)), ':g')
    % Starting time of Cue and CF
    time_cue = starting_time_cue(i)/512;
    time_cf = starting_time_cf(i)/512;
    xline(time_cue, 'k--')
    xline(time_cf, 'k--')
    name = char(subj_name);
    title(strcat('ERD avg logBP \beta - C4 - ',name(1:3)))
    legend('both feet', 'both hands')
    axis tight
    drawnow
    hold off

    % save average ERD in the structure
    subjects.(subj_name).ERD_logBP_mu_avg_feet    = ERD_logBP_mu_avg_feet;
    subjects.(subj_name).ERD_logBP_mu_avg_hands   = ERD_logBP_mu_avg_hands;
    subjects.(subj_name).ERD_logBP_beta_avg_feet  = ERD_logBP_beta_avg_feet;
    subjects.(subj_name).ERD_logBP_beta_avg_hands = ERD_logBP_beta_avg_hands;

end


%% Grand average for ERD on logarithmic band power

% Adjust different lengths
minLen_ERDmu_feet    = inf;
minLen_ERDmu_hands   = inf;
minLen_ERDbeta_feet  = inf;
minLen_ERDbeta_hands = inf;

for i = 1:length(data)
    subj_name = data(i);
    if length(subjects.(subj_name).ERD_logBP_mu_avg_feet) <  minLen_ERDmu_feet
        minLen_ERDmu_feet = length(subjects.(subj_name).ERD_logBP_mu_avg_feet);
    end
    if length(subjects.(subj_name).ERD_logBP_mu_avg_hands) <  minLen_ERDmu_hands
        minLen_ERDmu_hands = length(subjects.(subj_name).ERD_logBP_mu_avg_hands);
    end
    if length(subjects.(subj_name).ERD_logBP_beta_avg_feet) <  minLen_ERDbeta_feet
        minLen_ERDbeta_feet = length(subjects.(subj_name).ERD_logBP_mu_avg_feet);
    end
    if length(subjects.(subj_name).ERD_logBP_beta_avg_hands) <  minLen_ERDbeta_hands
        minLen_ERDbeta_hands = length(subjects.(subj_name).ERD_logBP_beta_avg_hands);
    end
end

% Save data from all subjects
ERDmu_feet_tot    = zeros(minLen_ERDmu_feet, 16, length(data));
ERDmu_hands_tot   = zeros(minLen_ERDmu_hands, 16, length(data));
ERDbeta_feet_tot  = zeros(minLen_ERDbeta_feet, 16, length(data));
ERDbeta_hands_tot = zeros(minLen_ERDbeta_hands, 16, length(data));

for i = 1:length(data)
    subj_name = data(i);
    ERDmu_feet_tot(:, :, i)    = subjects.(subj_name).ERD_logBP_mu_avg_feet(1:minLen_ERDmu_feet,:);
    ERDmu_hands_tot(:, :, i)   = subjects.(subj_name).ERD_logBP_mu_avg_hands(1:minLen_ERDmu_hands,:);
    ERDbeta_feet_tot(:, :, i)  = subjects.(subj_name).ERD_logBP_beta_avg_feet(1:minLen_ERDbeta_feet,:);
    ERDbeta_hands_tot(:, :, i) = subjects.(subj_name).ERD_logBP_beta_avg_hands(1:minLen_ERDbeta_hands,:);

end

% Compute Grand Average for all subjects
ERDmu_feet_GA    = mean(ERDmu_feet_tot, 3);
ERDmu_hands_GA   = mean(ERDmu_hands_tot, 3);
ERDbeta_feet_GA  = mean(ERDbeta_feet_tot, 3);
ERDbeta_hands_GA = mean(ERDbeta_hands_tot, 3);

% Compute standard error
ERDmu_feet_SE    = std(ERDmu_feet_tot, 0, 3)./sqrt(length(ERDmu_feet_tot(1, 1, :)));
ERDmu_hands_SE   = std(ERDmu_hands_tot, 0, 3)./sqrt(length(ERDmu_hands_tot(1, 1, :)));
ERDbeta_feet_SE  = std(ERDbeta_feet_tot, 0, 3)./sqrt(length(ERDbeta_feet_tot(1, 1, :)));
ERDbeta_hands_SE = std(ERDbeta_hands_tot, 0, 3)./sqrt(length(ERDbeta_hands_tot(1, 1, :)));

% GA temporal plots
hf = figure('units','normalized','outerposition',[0 0 1 1]); 
hf.Name = 'GA temporal plots';
hf.NumberTitle = 'off';
T = 1/sample_rate;                              % [s]
t = 0:T:(length(ERDmu_feet_GA(:,1,1))-1)*T;     % [s]

subplot(231), hold on
plot(t, ERDmu_feet_GA(:, chns(1)), 'g')
plot(t, ERDmu_feet_GA(:, chns(1)) - ERDmu_feet_SE(:, chns(1)), ':g')
plot(t, ERDmu_feet_GA(:, chns(1)) + ERDmu_feet_SE(:, chns(1)), ':g')
plot(t, ERDmu_hands_GA(:, chns(1)), 'r')
plot(t, ERDmu_hands_GA(:, chns(1)) - ERDmu_hands_SE(:, chns(1)), ':r')
plot(t, ERDmu_hands_GA(:, chns(1)) + ERDmu_hands_SE(:, chns(1)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
title('Grand Average ERD logBP \mu C3')
legend('both feet', 'both hands')
ylim([-40, 100])
hold off

subplot(232), hold on
plot(t, ERDmu_feet_GA(:, chns(2)), 'g')
plot(t, ERDmu_feet_GA(:, chns(2)) - ERDmu_feet_SE(:, chns(2)), ':g')
plot(t, ERDmu_feet_GA(:, chns(2)) + ERDmu_feet_SE(:, chns(2)), ':g')
plot(t, ERDmu_hands_GA(:, chns(2)), 'r')
plot(t, ERDmu_hands_GA(:, chns(2)) - ERDmu_hands_SE(:, chns(2)), ':r')
plot(t, ERDmu_hands_GA(:, chns(2)) + ERDmu_hands_SE(:, chns(2)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
ylim([-40, 100])
title('Grand Average ERD logBP \mu Cz')
legend('both feet', 'both hands')
hold off

subplot(233), hold on
plot(t, ERDmu_feet_GA(:, chns(3)), 'g')
plot(t, ERDmu_feet_GA(:, chns(3)) - ERDmu_feet_SE(:, chns(3)), ':g')
plot(t, ERDmu_feet_GA(:, chns(3)) + ERDmu_feet_SE(:, chns(3)), ':g')
plot(t, ERDmu_hands_GA(:, chns(3)), 'r')
plot(t, ERDmu_hands_GA(:, chns(3)) - ERDmu_hands_SE(:, chns(3)), ':r')
plot(t, ERDmu_hands_GA(:, chns(3)) + ERDmu_hands_SE(:, chns(3)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
ylim([-40, 100])
title('Grand Average ERD logBP \mu C4')
legend('both feet', 'both hands')
hold off

subplot(234), hold on
plot(t, ERDbeta_feet_GA(:, chns(1)), 'g')
plot(t, ERDbeta_feet_GA(:, chns(1)) - ERDbeta_feet_SE(:, chns(1)), ':g')
plot(t, ERDbeta_feet_GA(:, chns(1)) + ERDbeta_feet_SE(:, chns(1)), ':g')
plot(t, ERDbeta_hands_GA(:, chns(1)), 'r')
plot(t, ERDbeta_hands_GA(:, chns(1)) - ERDbeta_hands_SE(:, chns(1)), ':r')
plot(t, ERDbeta_hands_GA(:, chns(1)) + ERDbeta_hands_SE(:, chns(1)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
ylim([-40, 100])
title('Grand Average ERD logBP \beta C3')
legend('both feet', 'both hands')
hold off

subplot(235), hold on
plot(t, ERDbeta_feet_GA(:, chns(2)), 'g')
plot(t, ERDbeta_feet_GA(:, chns(2)) - ERDbeta_feet_SE(:, chns(2)), ':g')
plot(t, ERDbeta_feet_GA(:, chns(2)) + ERDbeta_feet_SE(:, chns(2)), ':g')
plot(t, ERDbeta_hands_GA(:, chns(2)), 'r')
plot(t, ERDbeta_hands_GA(:, chns(2)) - ERDbeta_hands_SE(:, chns(2)), ':r')
plot(t, ERDbeta_hands_GA(:, chns(2)) + ERDbeta_hands_SE(:, chns(2)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
ylim([-40, 100])
title('Grand Average ERD logBP \beta Cz')
legend('both feet', 'both hands')
hold off

subplot(236), hold on
plot(t, ERDbeta_feet_GA(:, chns(3)), 'g')
plot(t, ERDbeta_feet_GA(:, chns(3)) - ERDbeta_feet_SE(:, chns(3)), ':g')
plot(t, ERDbeta_feet_GA(:, chns(3)) + ERDbeta_feet_SE(:, chns(3)), ':g')
plot(t, ERDbeta_hands_GA(:, chns(3)), 'r')
plot(t, ERDbeta_hands_GA(:, chns(3)) - ERDbeta_hands_SE(:, chns(3)), ':r')
plot(t, ERDbeta_hands_GA(:, chns(3)) + ERDbeta_hands_SE(:, chns(3)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
ylim([-40, 100])
title('Grand Average ERD logBP \beta C4')
legend('both feet', 'both hands')
hold off


% Topographic plots
load('chanlocs16.mat');
len = min(h.DUR(h.TYP == 786));

% mu band
ERD_Ref_773 = mean(ERDmu_hands_GA(1:len, :), 1);
ERD_Act_773 = mean(ERDmu_hands_GA(len+1:end, :), 1);
ERD_Ref_771 = mean(ERDmu_feet_GA(1:len, :), 1);
ERD_Act_771 = mean(ERDmu_feet_GA(len+1:end, :), 1);

hf = figure('units','normalized','outerposition',[0 0 1 1]); 
hf.Name = 'Topographical maps mu';
hf.NumberTitle = 'off';
subplot(221)
topoplot(squeeze(ERD_Ref_773), chanlocs16);
title('Reference - \mu band - both hands')
colorbar
clim([-30, 60])
subplot(222)
topoplot(squeeze(ERD_Act_773), chanlocs16);
title('Activity - \mu band - both hands')
colorbar
clim([-30, 60])
subplot(223)
topoplot(squeeze(ERD_Ref_771), chanlocs16);
title('Reference - \mu band - both feet')
colorbar
clim([-30, 60])
subplot(224)
topoplot(squeeze(ERD_Act_771), chanlocs16);
title('Activity - \mu band - both feet')
colorbar
clim([-30, 60])

% beta band
ERD_Ref_773 = mean(ERDbeta_hands_GA(1:len, :), 1);
ERD_Act_773 = mean(ERDbeta_hands_GA(len+1:end, :), 1);
ERD_Ref_771 = mean(ERDbeta_feet_GA(1:len, :), 1);
ERD_Act_771 = mean(ERDbeta_feet_GA(len+1:end, :), 1);

hf = figure('units','normalized','outerposition',[0 0 1 1]); 
hf.Name = 'Topographical maps beta';
hf.NumberTitle = 'off';
subplot(221)
topoplot(squeeze(ERD_Ref_773), chanlocs16);
title('Reference - \beta band - both hands')
colorbar
clim([-15, 30])
subplot(222)
topoplot(squeeze(ERD_Act_773), chanlocs16);
title('Activity - \beta band - both hands')
colorbar
clim([-15, 30])
subplot(223)
topoplot(squeeze(ERD_Ref_771), chanlocs16);
title('Reference - \beta band - both feet')
colorbar
clim([-15, 30])
subplot(224)
topoplot(squeeze(ERD_Act_771), chanlocs16);
title('Activity - \beta band - both feet')
colorbar
clim([-15, 30])


%% ERD/ERS on spectrogram

% Concatenate the files

for i = 1:length(data)

    disp(['Subject ', data{i}(1:3)])

    subj_name = data(i);
    runs_names = fieldnames(subjects.(subj_name));

    % Initialize structure fields
    subjects.(subj_name).PSD_c = [];
    POS = [];
    DUR = [];

    % Concatenation of the PSD for the offline files
    for j = 1:length(runs_names)

        if contains(runs_names{j}, 'offline', 'IgnoreCase', true)

            DUR = [DUR; subjects.(subj_name).(runs_names{j}).h_PSD.DUR];
            POS = [POS; subjects.(subj_name).(runs_names{j}).h_PSD.POS + length(subjects.(subj_name).PSD_c)];

            subjects.(subj_name).PSD_c = [subjects.(subj_name).PSD_c; subjects.(subj_name).(runs_names{j}).PSD];        
        end

    end

    % Save information
    subjects.(subj_name).h_PSD.POS = POS;
    subjects.(subj_name).h_PSD.DUR = DUR;
    subjects.(subj_name).h_PSD.TYP = subjects.(subj_name).h.TYP;

    subjects.(subj_name).vectors_PSD = labelVecs(subjects.(subj_name).PSD_c, subjects.(subj_name).h_PSD);

    % Concatenation of the PSD for the online files
    subjects.(subj_name).PSD_c_online = [];
    TYP = [];
    POS = [];
    DUR = [];

    for j = 1:length(runs_names)

        if contains(runs_names{j}, 'online', 'IgnoreCase', true)

            TYP = [TYP; subjects.(subj_name).(runs_names{j}).h_PSD.TYP];
            DUR = [DUR; subjects.(subj_name).(runs_names{j}).h_PSD.DUR];
            POS = [POS; subjects.(subj_name).(runs_names{j}).h_PSD.POS + length(subjects.(subj_name).PSD_c_online)];

            subjects.(subj_name).PSD_c_online = [subjects.(subj_name).PSD_c_online; subjects.(subj_name).(runs_names{j}).PSD];        
        end

    end

    subjects.(subj_name).h_PSD_online.POS = POS;
    subjects.(subj_name).h_PSD_online.DUR = DUR;
    subjects.(subj_name).h_PSD_online.TYP = TYP;

    subjects.(subj_name).vectors_PSD_online = labelVecs(subjects.(subj_name).PSD_c_online, subjects.(subj_name).h_PSD_online);
    
end

clc
disp('Done')

%% Activity and Reference computation

ax1 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax1.Name = 'Average ERD for both feet task - channel C3'; % 7
ax1.NumberTitle = 'off';
hold off

ax2 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax2.Name = 'Average ERD for both feet task - channel Cz'; % 9
ax2.NumberTitle = 'off';
hold off

ax3 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax3.Name = 'Average ERD for both feet task - channel C4'; % 11
ax3.NumberTitle = 'off';
hold off

ax4 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax4.Name = 'Average ERD for both hands task - channel C3'; % 7
ax4.NumberTitle = 'off';
hold off

ax5 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax5.Name = 'Average ERD for both hands task - channel Cz'; % 9
ax5.NumberTitle = 'off';
hold off

ax6 = figure('units','normalized','outerposition',[0 0 1 1]); 
ax6.Name = 'Average ERD for both hands task - channel C4'; % 11
ax6.NumberTitle = 'off';
hold off

for i = 1:length(data)
    subj_name = data(i);
    runs_names = fieldnames(subjects.(subj_name));

    subjects.(subj_name).ERD = [];

    % Store the values in temporary variables
    h_PSD = subjects.(subj_name).h_PSD;
    PSD_c = subjects.(subj_name).PSD_c;

    % Get the starting and ending positions of the trials
    startTrial = h_PSD.POS(h_PSD.TYP == 786);
    stopTrial  = h_PSD.POS(h_PSD.TYP == 781) + h_PSD.DUR(h_PSD.TYP == 781);

    % Get the number and length of trials
    ntrials = length(startTrial);
    trial_length = min(stopTrial - startTrial);

    % Initialize data strucure
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

    ERD = 100 * (Activity- Reference)./ Reference;
    % ERD = log(Activity ./ Reference);
    % ERD = 10 * log10(Activity ./ Reference);
    
    subjects.(subj_name).ERD_PSD = ERD;

    ERDavg_feet  = mean(ERD(:, :, :, subjects.(subj_name).vectors_PSD.Ck == 771), 4);
    ERDavg_hands = mean(ERD(:, :, :, subjects.(subj_name).vectors_PSD.Ck == 773), 4);

    subjects.(subj_name).ERDavg_feet  = ERDavg_feet;
    subjects.(subj_name).ERDavg_hands = ERDavg_hands;

    % Visualization
    
    % frequency vector
    f = subjects.(subj_name).(runs_names{1}).f;     % [Hz]
    % time vector
    T = 0.0627;                                     % wshift, [s]
    t = 0:T:length(ERDavg_hands)*T;                 % [s]

    set(0,'CurrentFigure',ax1)
    subplot(2, 4, mod(i-1, 8)+1)
    imagesc(t, f, ERDavg_feet(:, :, 7)')
    colormap hot
    colorbar
    name = char(subj_name);
    title(strcat('ERD avg both feet - C3 - ',name(1:3)))
    set(gca,'YDir','normal')
    drawnow
    hold off

    set(0,'CurrentFigure',ax2)
    subplot(2, 4, mod(i-1, 8)+1, 'Parent', ax2)
    imagesc(t, f, ERDavg_feet(:, :, 9)')
    colormap hot
    colorbar
    name = char(subj_name);
    title(strcat('ERD avg both feet - Cz - ',name(1:3)))
    set(gca,'YDir','normal')
    drawnow
    hold off

    set(0,'CurrentFigure',ax3)
    subplot(2, 4, mod(i-1, 8)+1, 'Parent', ax3)
    imagesc(t, f, ERDavg_feet(:, :, 11)')
    colormap hot
    colorbar
    name = char(subj_name);
    title(strcat('ERD avg both feet - C4 - ',name(1:3)))
    set(gca,'YDir','normal')
    drawnow
    hold off

    set(0,'CurrentFigure',ax4)
    subplot(2, 4, mod(i-1, 8)+1, 'Parent', ax4)
    imagesc(t, f, ERDavg_hands(:, :, chns(1))')
    colormap hot
    colorbar
    name = char(subj_name);
    title(strcat('ERD avg both hands - C3 - ',name(1:3)))
    set(gca,'YDir','normal')
    drawnow
    hold off

    set(0,'CurrentFigure',ax5)
    subplot(2, 4, mod(i-1, 8)+1, 'Parent', ax5)
    imagesc(t, f, ERDavg_hands(:, :, chns(2))')
    colormap hot
    colorbar
    name = char(subj_name);
    title(strcat('ERD avg both hands - Cz - ',name(1:3)))
    set(gca,'YDir','normal')
    drawnow
    hold off

    set(0,'CurrentFigure',ax6)
    subplot(2, 4, mod(i-1, 8)+1, 'Parent', ax6)
    imagesc(t, f, ERDavg_hands(:, :, chns(3))')
    colormap hot
    colorbar
    name = char(subj_name);
    title(strcat('ERD avg both hands - C4 - ',name(1:3)))
    set(gca,'YDir','normal')
    drawnow
    hold off

end

%% Grand average ERD/ERS on PSD of all subjects

% Adjust different lengths
minLen_ERDfeet    = inf;
minLen_ERDhands   = inf;

for i = 1:length(data)
    subj_name = data(i);
    if length(subjects.(subj_name).ERDavg_feet) <  minLen_ERDfeet
        minLen_ERDfeet = length(subjects.(subj_name).ERDavg_feet);
    end
    if length(subjects.(subj_name).ERDavg_hands) <  minLen_ERDhands
        minLen_ERDhands = length(subjects.(subj_name).ERDavg_hands);
    end
end

% Initialize data structures
ERDfeet_tot  = zeros(minLen_ERDfeet, length(f), length(channels), length(data));
ERDhands_tot = zeros(minLen_ERDhands, length(f), length(channels), length(data));

for i = 1:length(data)
    subj_name = data(i);
    ERDfeet_tot(:, :, :, i)    = subjects.(subj_name).ERDavg_feet(1:minLen_ERDfeet,:,:);
    ERDhands_tot(:, :, :, i)   = subjects.(subj_name).ERDavg_hands(1:minLen_ERDhands,:,:);
end

% Compute Grand Average for all subjects
ERDfeet_GA  = mean(ERDfeet_tot, 4);
ERDhands_GA = mean(ERDhands_tot, 4);

hf = figure('units','normalized','outerposition',[0 0 1 1]); 
hf.Name = 'GA ERD on PSD';
hf.NumberTitle = 'off';
subplot(231)
imagesc(t, f, ERDhands_GA(:, :, chns(1))')
colormap hot
colorbar
title('GA ERD avg both hands - C3')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(232)
imagesc(t, f, ERDhands_GA(:, :, chns(2))')
colormap hot
colorbar
title('GA ERD avg both hands - Cz')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(233)
imagesc(t, f, ERDhands_GA(:, :, chns(3))')
colormap hot
colorbar
title('GA ERD avg both hands - C4')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(234)
imagesc(t, f, ERDfeet_GA(:, :, chns(1))')
colormap hot
colorbar
title('GA ERD avg both feet - C3')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(235)
imagesc(t, f, ERDfeet_GA(:, :, chns(2))')
colormap hot
colorbar
title('GA ERD avg both feet - Cz')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(236)
imagesc(t, f, ERDfeet_GA(:, :, chns(3))')
colormap hot
colorbar
title('GA ERD avg both feet - C4')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')


%% Features calculation and maps visualization on each run

for i = 1:length(data)
    subj_name = data(i);

    runs = dir(fullfile(pwd, strcat("Data/", subj_name)));
    runs_names = {runs.name};
    
    count_off = 0;
    for j = 1:length(runs_names)
        run_name = runs_names{j};
        if run_name(1) == 'a'

            % Offline runs
            if run_name(21:27) == 'offline'
                count_off = count_off + 1;
                field_name = strcat("offline",string(count_off));

                % Store the values in temporary variables
                PSD = subjects.(subj_name).(field_name).PSD;
                h_PSD = subjects.(subj_name).(field_name).h_PSD;
                h_PSD.TYP = subjects.(subj_name).(field_name).h.EVENT.TYP;      % Add TYP field

                % Calculate the label vectors
                tmpVecs = labelVecs(PSD, h_PSD);

                % Extract PSD windows from Cue to Continuous Feedback and take the logarithm
                wnds_CktoCFk = log(PSD(tmpVecs.Ak > 0 | tmpVecs.CFk > 0, :, :));

                % Reshape features to obtain a 2D matrix [windows x (channels x frequencies)]
                features = reshape(wnds_CktoCFk, [size(wnds_CktoCFk, 1), size(wnds_CktoCFk, 2) * size(wnds_CktoCFk, 3)]);

                idx = tmpVecs.Ak + tmpVecs.CFk;                 % Vector containing the codes for Cue and Continuous Feedback
                idx = idx(idx > 0);
                class = (idx < 781) .* idx;                     % Vector containing the class for each window

                for cl = 2 : length(idx)
                    if class(cl-1) > 0 && class(cl) == 0        % If the current class is 0, then
                        class(cl) = class(cl-1);                % the value is updated according
                    end                                         % to the previous state.
                end

                % Fisher score computation
                FS = abs(mean(features(class == 771, :), 1) - mean(features(class == 773, :), 1))./sqrt(std(features(class == 771, :), 1).^2 + std(features(class == 773, :), 1).^2);

                % Feature maps visualization
                hf = figure;
                hf.Name = ['Subject ', data{i}(1:3), ' run ', num2str(j-2)];
                hf.NumberTitle = 'off';
                title(['Subject ', data{i}(1:3), ' run ', num2str(j-2)])
                imagesc(f, 1:16, flipud(imrotate(reshape(FS, [23, 16]), 90)))       % Visualize Fisher Score maps
                xticks(f)
                xtickangle(90)
                xlabel("Hz")
                set(gca,'YTick', 1:16)
                set(gca,'yticklabel',(channels))
                ylabel("Channel")
                hold on

                [row_feat, col_feat] = find(ismember(reshape(FS, [23, 16]), maxk(FS, 3)));      % Indeces of the 3 most discriminative features

                plot(f(row_feat), col_feat, 'ro', 'MarkerSize', 10, 'LineWidth', 2)             % Circle the most discriminative features
                hold off

            end
        end
    end
end

%% Features calculation and maps visualization on concatenated runs


ax0 = figure('units','normalized','outerposition',[0 0 1 1]); 

for i = 1:length(data)
    subj_name = data(i);

    % Extract PSD windows from Cue to Continuous Feedback and take the logarithm
    wnds_CktoCFk = log(subjects.(subj_name).PSD_c(subjects.(subj_name).vectors_PSD.Ak > 0 | subjects.(subj_name).vectors_PSD.CFk > 0, :, :));
    
    % Reshape features to obtain a 2D matrix [windows x (channels x frequencies)]
    features = reshape(wnds_CktoCFk, [size(wnds_CktoCFk, 1), size(wnds_CktoCFk, 2) * size(wnds_CktoCFk, 3)]);

    idx = subjects.(subj_name).vectors_PSD.Ak + subjects.(subj_name).vectors_PSD.CFk;       % Vector containing the codes for Cue and Continuous Feedback
    idx = idx(idx > 0);
    class = (idx < 781) .* idx;                                                             % Vector containing the class for each window

    for j = 2 : length(idx)
        if class(j-1) > 0 && class(j) == 0      % If the current class is 0, then
            class(j) = class(j-1);              % the value is updated according
        end                                     % to the previous state.
    end

    % Fisher score computation
    FS = abs(mean(features(class == 771, :), 1) - mean(features(class == 773, :), 1))./sqrt(std(features(class == 771, :), 1).^2 + std(features(class == 773, :), 1).^2);

    % Feature maps visualization
    set(0,'CurrentFigure',ax0)
    subplot(2, 4, mod(i-1, 8)+1, 'Parent', ax0)
    title(strcat('Subject ', data{i}(1:3)))
    imagesc(f, 1:16, flipud(imrotate(reshape(FS, [23, 16]), 90)))
    xticks(f)
    xtickangle(90)
    xlabel("Hz")
    set(gca,'YTick', 1:16)
    set(gca,'yticklabel',(channels))
    ylabel("Channel")
    hold on

    [row_feat, col_feat] = find(ismember(reshape(FS, [23, 16]), maxk(FS, 3)));      % Indeces of the 3 most discriminative features

    % Save important information for the models (offline runs)
    subjects.(subj_name).wnds_CktoCFk = wnds_CktoCFk;
    subjects.(subj_name).class = class;
    subjects.(subj_name).row_feat = row_feat;
    subjects.(subj_name).col_feat = col_feat;

    plot(f(row_feat), col_feat, 'ro', 'MarkerSize', 10, 'LineWidth', 2)             % Circle the most discriminative features
    hold off

    % Define features for online runs
    wnds_CktoCFk_online = log(subjects.(subj_name).PSD_c_online(subjects.(subj_name).vectors_PSD_online.Ak > 0 | subjects.(subj_name).vectors_PSD_online.CFk > 0, :, :));
    features_online = reshape(wnds_CktoCFk_online, [size(wnds_CktoCFk_online, 1), size(wnds_CktoCFk_online, 2) * size(wnds_CktoCFk_online, 3)]);

    idx = subjects.(subj_name).vectors_PSD_online.Ak + subjects.(subj_name).vectors_PSD_online.CFk;
    idx = idx(idx > 0);
    class = (idx ~= 781) .* idx;

    for j = 2 : length(idx)
        if class(j-1) > 0 && class(j) == 0      % If the current class is 0, then
            class(j) = class(j-1);              % the value is updated according
        end                                     % to the previous state.
    end

    % Save important information (online runs)
    subjects.(subj_name).wnds_CktoCFk_online = wnds_CktoCFk_online;
    subjects.(subj_name).features_online = features_online;
    subjects.(subj_name).class_online = class;

end


%% Create classifiers based on the extracted features

for i = 1:length(data)
    subj_name = data(i);

    % Define the new matrix containing the data corresponding to the most discriminative features
    row_feat = subjects.(subj_name).row_feat;
    col_feat = subjects.(subj_name).col_feat;

    train_set = zeros(size(subjects.(subj_name).wnds_CktoCFk, 1), length(col_feat));

    for j = 1 : length(col_feat)
        train_set(:, j) = subjects.(subj_name).wnds_CktoCFk(:, row_feat(j), col_feat(j));
    end

    subjects.(subj_name).train_set = train_set;

    X = train_set;
    y = subjects.(subj_name).class;

    fprintf('Training the model for subject %s\n', data{i}(1:3))

    mdl = fitcdiscr(X, y, 'DiscrimType','quadratic');
    subjects.(subj_name).mdl = mdl;

    % Store all important information for each subject
    data_subjects.PSD_c_online = subjects.(subj_name).PSD_c_online;
    data_subjects.h_PSD_online = subjects.(subj_name).h_PSD_online;
    data_subjects.vectors_PSD_online = subjects.(subj_name).vectors_PSD_online;
    data_subjects.wnds_CktoCFk_online = subjects.(subj_name).wnds_CktoCFk_online;
    data_subjects.features_online = subjects.(subj_name).features_online;
    data_subjects.row_feat = subjects.(subj_name).row_feat;
    data_subjects.col_feat = subjects.(subj_name).col_feat;
    data_subjects.class = subjects.(subj_name).class_online;

    name = strcat(pwd, '\Data\', subj_name, '\data');
    
    % Save information in a .mat file
    save(name, 'mdl', 'data_subjects')
    
end

%% Models evaluation on training data

for i = 1:length(data)
    subj_name = data(i);

    mdl = subjects.(subj_name).mdl;

    % Model prediction
    [Gk, pp] = predict(mdl, subjects.(subj_name).train_set);

    % Save raw probability
    subjects.(subj_name).pp = pp;

    % Load the label vector
    y = subjects.(subj_name).class;

    % Calculate the total accuracy (correct predictions over total predictions)
    tot_accuracy = mean(Gk == y) * 100;
    
    % Calculate class accuracy separately
    feet_accuracy  = mean(Gk(y == 771) == y(y == 771)) * 100;
    hands_accuracy = mean(Gk(y == 773) == y(y == 773)) * 100;
    
    % Print accuracies
    fprintf('Accuracies of the model for subject %s\n', data{i}(1:3));
    fprintf('Accuracy: %f%%\n', tot_accuracy);
    fprintf('Accuracy both feet: %f%%\n', feet_accuracy);
    fprintf('Accuracy both hands: %f%%\n\n', hands_accuracy);

    % Accuracy bar plot
    xaxis = ["overall" "both hands" "both feet"];
    accuracies = [tot_accuracy, hands_accuracy, feet_accuracy];

    % Bar graphs
    figure
    title(['Single sample accuracy on training set for subject ', data{i}(i)])
    bar(xaxis, accuracies)
    ylabel('Accuracy [%]')
    ylim([0, 100])
    grid on
end

%% Accumulation framework for offline data

% Exponential accumulation framework
alpha = 0.95;                   % Smoothing parameter

for i = 1:length(data)
    subj_name = data(i);
    pp = subjects.(subj_name).pp;

    % For each window in Cue and Continuous Feedback tasks, extract the number of the trial
    trials_windows = subjects.(subj_name).vectors_PSD.Tk(subjects.(subj_name).vectors_PSD.Ak > 0 | subjects.(subj_name).vectors_PSD.CFk > 0);
    
    % Define the number of windows
    nwindows = length(pp);

    % Initialize Decision data structure
    D = 0.5 * ones(nwindows, 2);

    for wId = 2:nwindows

        % Verify if the sample is the start of a new trial
        if trials_windows(wId) ~= trials_windows(wId-1)
            % If it is, restore the Decision value to 0.5
            D(wId,:) = [0.5 0.5];
        else 
            % Otherwise, compute the Decision value based on pp and previous Decision
            D(wId,:) = D(wId - 1,:) * alpha + pp(wId,:) * (1 - alpha);
        end
    end
    
    % Save results
    subjects.(subj_name).D = D;

end

%% Trial Accuracy computation

% Define decision thresholds
thr = [0.3 0.7];

% Initialize data structure to store the average time to deliver a command 
t_avg = zeros(1, length(data));

for i = 1:length(data)
    subj_name = data(i);

    % For each window in Cue and Continuous Feedback tasks, extract the number of the trial
    trials_windows = subjects.(subj_name).vectors_PSD.Tk(subjects.(subj_name).vectors_PSD.Ak > 0 | subjects.(subj_name).vectors_PSD.CFk > 0);
    pp = subjects.(subj_name).pp;
    D = subjects.(subj_name).D;

    % Initialize structures
    Gk_trial_all = zeros(max(trials_windows),1);            % Store trial prediction for each trial
    t_for_command = zeros(max(trials_windows),1);           % Store the time for a command for each trial

    % Compute classification for succesful trials
    for trial_number = 1 : max(trials_windows)
        
        D_trial = D(trials_windows == trial_number,:);

        % Structure to store the sample when the threshold is reached
        t = zeros(1, length(D_trial));

        for j = 1 : length(D_trial)
            if D_trial(j) <= thr(1)
                Gk_trial_all(trial_number) = 773;   % Trial classified as both hands
                t(j) = j;                           % Save the sample index if the th is reached
            end
            if D_trial(j) >= thr(2)
                Gk_trial_all(trial_number) = 771;   % Trial classified as both feet
                t(j) = j;                           % Save the sample index if the th is reached

            end
        end

        % If the threshold is reached, t_for_command is saved as the first sample after reaching the threshold
        if sum(t) > 0
            t_for_command(trial_number) = find(t > 0, 1, "first");
        end
    end

    % Compute the average time to deliver a command
    t_avg(i) = mean(t_for_command(t_for_command > 0));

    % Save trial predictions
    subjects.(subj_name).Gk_trial_all = Gk_trial_all;
end

% Convert time in seconds
t_avg = t_avg * 0.0625;             % [s]


%% Trial accuracy with and without rejection of time-out trials

% Clear command window to print new accuracies
clc

% Initialize data structures
trial_accuracy_no_rejection = zeros(1, length(data));
trial_accuracy_no_rejection_feet  = zeros(1, length(data));
trial_accuracy_no_rejection_hands = zeros(1, length(data));

for i = 1:length(data)
    subj_name = data(i);

    % Store values in temporary variables
    trials_windows = subjects.(subj_name).vectors_PSD.Tk(subjects.(subj_name).vectors_PSD.Ak > 0 | subjects.(subj_name).vectors_PSD.CFk > 0);
    pp = subjects.(subj_name).pp;
    D = subjects.(subj_name).D;
    Gk_trial_all = subjects.(subj_name).Gk_trial_all;

    % Calculate accuracy
    trial_accuracy_no_rejection(i) = mean(Gk_trial_all == subjects.(subj_name).vectors_PSD.Ck) * 100;
    trial_accuracy_no_rejection_feet(i)  = mean(Gk_trial_all(subjects.(subj_name).vectors_PSD.Ck == 771) == subjects.(subj_name).vectors_PSD.Ck(subjects.(subj_name).vectors_PSD.Ck == 771)) * 100;
    trial_accuracy_no_rejection_hands(i) = mean(Gk_trial_all(subjects.(subj_name).vectors_PSD.Ck == 773) == subjects.(subj_name).vectors_PSD.Ck(subjects.(subj_name).vectors_PSD.Ck == 773)) * 100;

    % Data for the bar graph
    accuracies = [trial_accuracy_no_rejection(i), trial_accuracy_no_rejection_feet(i), trial_accuracy_no_rejection_hands(i)];
    x_labels = {'Overall', 'Both Hands', 'Both Feet'};
    
    % Accuracies visualization as a bar graph
    hf = figure;
    hf.Name = ['Subject ', data{i}(1:3), ' accuracies on test'];
    hf.NumberTitle = 'off';
    bar(accuracies)
    set(gca, 'xticklabel', x_labels)
    ylabel('Accuracy [%]')
    title(['Trial accuracy on test set - subject ', data{i}(1:3)])
    ylim([0, 100])
    grid on

    % Print accuracies
    fprintf('Trial accuracies of the model (no rejection) for subject %s\n', data{i}(1:3));
    fprintf('Accuracy: %f%%\n', trial_accuracy_no_rejection(i));
    fprintf('Accuracy both feet: %f%%\n', trial_accuracy_no_rejection_feet(i));
    fprintf('Accuracy both hands: %f%%\n\n', trial_accuracy_no_rejection_hands(i));
end

%% Rejection

% Clear command window to print new accuracies
clc

% Initialize data structures
trial_accuracy_rejection = zeros(1, length(data));
trial_accuracy_rejection_feet  = zeros(1, length(data));
trial_accuracy_rejection_hands = zeros(1, length(data));

for i = 1:length(data)
    subj_name = data(i);

    % Store values in temporary variables
    trials_windows = subjects.(subj_name).vectors_PSD.Tk(subjects.(subj_name).vectors_PSD.Ak > 0 | subjects.(subj_name).vectors_PSD.CFk > 0);
    pp = subjects.(subj_name).pp;
    D = subjects.(subj_name).D;
    Gk_trial_all = subjects.(subj_name).Gk_trial_all;

    % Find the trials to reject
    idx_rejection = (Gk_trial_all ~= 0);

    % Calculate the accuracy on the remaining trials
    trial_accuracy_rejection(i) = mean(Gk_trial_all(idx_rejection) == subjects.(subj_name).vectors_PSD.Ck(idx_rejection)) * 100;
    trial_accuracy_rejection_feet(i)  = mean(Gk_trial_all(idx_rejection & subjects.(subj_name).vectors_PSD.Ck == 771) == subjects.(subj_name).vectors_PSD.Ck(idx_rejection & subjects.(subj_name).vectors_PSD.Ck == 771)) * 100;
    trial_accuracy_rejection_hands(i) = mean(Gk_trial_all(idx_rejection & subjects.(subj_name).vectors_PSD.Ck == 773) == subjects.(subj_name).vectors_PSD.Ck(idx_rejection & subjects.(subj_name).vectors_PSD.Ck == 773)) * 100;

    % Data for the bar graph
    accuracies = [trial_accuracy_rejection(i), trial_accuracy_rejection_feet(i), trial_accuracy_rejection_hands(i)];
    x_labels = {'Overall', 'Both Hands', 'Both Feet'};
    
    % Accuracies visualization as a bar graph
    hf = figure;
    hf.Name = ['Subject ', data{i}(1:3), ' accuracies on test'];
    hf.NumberTitle = 'off';
    bar(accuracies)
    set(gca, 'xticklabel', x_labels)
    ylabel('Accuracy [%]')
    title(['Trial accuracy on test set - subject ', data{i}(1:3)])
    ylim([0, 100])
    grid on

    % Print accuracies
    fprintf('Trial accuracies of the model for subject %s\n', data{i}(1:3));
    fprintf('Accuracy: %f%%\n', trial_accuracy_rejection(i));
    fprintf('Accuracy both feet: %f%%\n', trial_accuracy_rejection_feet(i));
    fprintf('Accuracy both hands: %f%%\n', trial_accuracy_rejection_hands(i));

    % Print time to deliver a command
    fprintf('Time to deliver a command %f s\n\n', t_avg(i));

end



%% Grand Average analysis on representative subjects

% In order to find representative subjects, each subject's feature maps
% were considered. A subject was considered as representative if:
% 1) The features of each of his runs were stable;
% 2) The features included (at least 2) of the most relevant channels for
% the analysis (C3, Cz, C4).

% For these reasons, subjects aj7, aj9 and ai7 were not considered as
% representative:
% Aj7: features instability during the runs
% Aj9: no relevant channel for the analysis among the features
% Ai7: features instability during the runs


% Subjects to remove based on the feature maps
subjects_to_remove = {'aj7', 'aj9', 'ai7'};

% Initialize data structures
ERDmu_feet_tot    = zeros(minLen_ERDmu_feet, length(channels), length(data)-length(subjects_to_remove));
ERDmu_hands_tot   = zeros(minLen_ERDmu_hands, length(channels), length(data)-length(subjects_to_remove));
ERDbeta_feet_tot  = zeros(minLen_ERDbeta_feet, length(channels), length(data)-length(subjects_to_remove));
ERDbeta_hands_tot = zeros(minLen_ERDbeta_hands, length(channels), length(data)-length(subjects_to_remove));

for i = 1 : (length(data) - length(subjects_to_remove))
    subj_name = data(i);

    % Check that the name is not associated with a non representative subject
    if ~contains(subj_name, 'aj7', 'IgnoreCase', true) && ~contains(subj_name, 'aj9', 'IgnoreCase', true) && ~contains(subj_name, 'ai7', 'IgnoreCase', true)
    % if ~contains(subj_name, 'aj7', 'IgnoreCase', true) && ~contains(subj_name, 'ai7', 'IgnoreCase', true)
        ERDmu_feet_tot(:, :, i)    = subjects.(subj_name).ERD_logBP_mu_avg_feet(1:minLen_ERDmu_feet,:);
        ERDmu_hands_tot(:, :, i)   = subjects.(subj_name).ERD_logBP_mu_avg_hands(1:minLen_ERDmu_hands,:);
        ERDbeta_feet_tot(:, :, i)  = subjects.(subj_name).ERD_logBP_beta_avg_feet(1:minLen_ERDbeta_feet,:);
        ERDbeta_hands_tot(:, :, i) = subjects.(subj_name).ERD_logBP_beta_avg_hands(1:minLen_ERDbeta_hands,:);
    end
end

% Compute Grand Average
ERDmu_feet_GA    = mean(ERDmu_feet_tot, 3);
ERDmu_hands_GA   = mean(ERDmu_hands_tot, 3);
ERDbeta_feet_GA  = mean(ERDbeta_feet_tot, 3);
ERDbeta_hands_GA = mean(ERDbeta_hands_tot, 3);

% Compute standard error
ERDmu_feet_SE    = std(ERDmu_feet_tot, 0, 3)./sqrt(length(ERDmu_feet_tot(1, 1, :)));
ERDmu_hands_SE   = std(ERDmu_hands_tot, 0, 3)./sqrt(length(ERDmu_hands_tot(1, 1, :)));
ERDbeta_feet_SE  = std(ERDbeta_feet_tot, 0, 3)./sqrt(length(ERDbeta_feet_tot(1, 1, :)));
ERDbeta_hands_SE = std(ERDbeta_hands_tot, 0, 3)./sqrt(length(ERDbeta_hands_tot(1, 1, :)));

% GA temporal plots
hf = figure('units','normalized','outerposition',[0 0 1 1]); 
hf.Name = 'GA temporal plots - representative subjects';
hf.NumberTitle = 'off';
T = 1/sample_rate;                              % [s]
t = 0:T:(length(ERDmu_feet_GA(:,1,1))-1)*T;     % [s]

subplot(231), hold on
plot(t, ERDmu_feet_GA(:, chns(1)), 'g')
plot(t, ERDmu_feet_GA(:, chns(1)) - ERDmu_feet_SE(:, chns(1)), ':g')
plot(t, ERDmu_feet_GA(:, chns(1)) + ERDmu_feet_SE(:, chns(1)), ':g')
plot(t, ERDmu_hands_GA(:, chns(1)), 'r')
plot(t, ERDmu_hands_GA(:, chns(1)) - ERDmu_hands_SE(:, chns(1)), ':r')
plot(t, ERDmu_hands_GA(:, chns(1)) + ERDmu_hands_SE(:, chns(1)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
title('Grand Average ERD logBP \mu C3')
legend('both feet', 'both hands')
ylim([-50, 70])
hold off

subplot(232), hold on
plot(t, ERDmu_feet_GA(:, chns(2)), 'g')
plot(t, ERDmu_feet_GA(:, chns(2)) - ERDmu_feet_SE(:, chns(2)), ':g')
plot(t, ERDmu_feet_GA(:, chns(2)) + ERDmu_feet_SE(:, chns(2)), ':g')
plot(t, ERDmu_hands_GA(:, chns(2)), 'r')
plot(t, ERDmu_hands_GA(:, chns(2)) - ERDmu_hands_SE(:, chns(2)), ':r')
plot(t, ERDmu_hands_GA(:, chns(2)) + ERDmu_hands_SE(:, chns(2)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
ylim([-50, 70])
title('Grand Average ERD logBP \mu Cz')
legend('both feet', 'both hands')
hold off

subplot(233), hold on
plot(t, ERDmu_feet_GA(:, chns(3)), 'g')
plot(t, ERDmu_feet_GA(:, chns(3)) - ERDmu_feet_SE(:, chns(3)), ':g')
plot(t, ERDmu_feet_GA(:, chns(3)) + ERDmu_feet_SE(:, chns(3)), ':g')
plot(t, ERDmu_hands_GA(:, chns(3)), 'r')
plot(t, ERDmu_hands_GA(:, chns(3)) - ERDmu_hands_SE(:, chns(3)), ':r')
plot(t, ERDmu_hands_GA(:, chns(3)) + ERDmu_hands_SE(:, chns(3)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
ylim([-50, 70])
title('Grand Average ERD logBP \mu C4')
legend('both feet', 'both hands')
hold off

subplot(234), hold on
plot(t, ERDbeta_feet_GA(:, chns(1)), 'g')
plot(t, ERDbeta_feet_GA(:, chns(1)) - ERDbeta_feet_SE(:, chns(1)), ':g')
plot(t, ERDbeta_feet_GA(:, chns(1)) + ERDbeta_feet_SE(:, chns(1)), ':g')
plot(t, ERDbeta_hands_GA(:, chns(1)), 'r')
plot(t, ERDbeta_hands_GA(:, chns(1)) - ERDbeta_hands_SE(:, chns(1)), ':r')
plot(t, ERDbeta_hands_GA(:, chns(1)) + ERDbeta_hands_SE(:, chns(1)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
ylim([-50, 70])
title('Grand Average ERD logBP \beta C3')
legend('both feet', 'both hands')
hold off

subplot(235), hold on
plot(t, ERDbeta_feet_GA(:, chns(2)), 'g')
plot(t, ERDbeta_feet_GA(:, chns(2)) - ERDbeta_feet_SE(:, chns(2)), ':g')
plot(t, ERDbeta_feet_GA(:, chns(2)) + ERDbeta_feet_SE(:, chns(2)), ':g')
plot(t, ERDbeta_hands_GA(:, chns(2)), 'r')
plot(t, ERDbeta_hands_GA(:, chns(2)) - ERDbeta_hands_SE(:, chns(2)), ':r')
plot(t, ERDbeta_hands_GA(:, chns(2)) + ERDbeta_hands_SE(:, chns(2)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
ylim([-50, 70])
title('Grand Average ERD logBP \beta Cz')
legend('both feet', 'both hands')
hold off

subplot(236), hold on
plot(t, ERDbeta_feet_GA(:, chns(3)), 'g')
plot(t, ERDbeta_feet_GA(:, chns(3)) - ERDbeta_feet_SE(:, chns(3)), ':g')
plot(t, ERDbeta_feet_GA(:, chns(3)) + ERDbeta_feet_SE(:, chns(3)), ':g')
plot(t, ERDbeta_hands_GA(:, chns(3)), 'r')
plot(t, ERDbeta_hands_GA(:, chns(3)) - ERDbeta_hands_SE(:, chns(3)), ':r')
plot(t, ERDbeta_hands_GA(:, chns(3)) + ERDbeta_hands_SE(:, chns(3)), ':r')
xlabel('Time [s]')
ylabel('[ERD/ERS]')
ylim([-50, 70])
title('Grand Average ERD logBP \beta C4')
legend('both feet', 'both hands')
hold off


% Topographic plots
load('chanlocs16.mat');
len = min(h.DUR(h.TYP == 786));

% mu band
ERD_Ref_773 = mean(ERDmu_hands_GA(1:len, :), 1);
ERD_Act_773 = mean(ERDmu_hands_GA(len+1:end, :), 1);
ERD_Ref_771 = mean(ERDmu_feet_GA(1:len, :), 1);
ERD_Act_771 = mean(ERDmu_feet_GA(len+1:end, :), 1);

hf = figure('units','normalized','outerposition',[0 0 1 1]); 
hf.Name = 'Topographical maps mu';
hf.NumberTitle = 'off';
subplot(221)
topoplot(squeeze(ERD_Ref_773), chanlocs16);
title('Reference - \mu band - both hands')
colorbar
clim([-20, 20])
subplot(222)
topoplot(squeeze(ERD_Act_773), chanlocs16);
title('Activity - \mu band - both hands')
colorbar
clim([-20, 20])
subplot(223)
topoplot(squeeze(ERD_Ref_771), chanlocs16);
title('Reference - \mu band - both feet')
colorbar
clim([-20, 20])
subplot(224)
topoplot(squeeze(ERD_Act_771), chanlocs16);
title('Activity - \mu band - both feet')
colorbar
clim([-20, 20])

% beta band
ERD_Ref_773 = mean(ERDbeta_hands_GA(1:len, :), 1);
ERD_Act_773 = mean(ERDbeta_hands_GA(len+1:end, :), 1);
ERD_Ref_771 = mean(ERDbeta_feet_GA(1:len, :), 1);
ERD_Act_771 = mean(ERDbeta_feet_GA(len+1:end, :), 1);

hf = figure('units','normalized','outerposition',[0 0 1 1]); 
hf.Name = 'Topographical maps beta';
hf.NumberTitle = 'off';
subplot(221)
topoplot(squeeze(ERD_Ref_773), chanlocs16);
title('Reference - \beta band - both hands')
colorbar
clim([-10, 10])
subplot(222)
topoplot(squeeze(ERD_Act_773), chanlocs16);
title('Activity - \beta band - both hands')
colorbar
clim([-10, 10])
subplot(223)
topoplot(squeeze(ERD_Ref_771), chanlocs16);
title('Reference - \beta band - both feet')
colorbar
clim([-10, 10])
subplot(224)
topoplot(squeeze(ERD_Act_771), chanlocs16);
title('Activity - \beta band - both feet')
colorbar
clim([-10, 10])


%% Grand average ERD/ERS on PSD of representative subjects

% Initialize data structures
ERDfeet_tot  = zeros(minLen_ERDfeet, length(f), length(channels), length(data) - length(subjects_to_remove));
ERDhands_tot = zeros(minLen_ERDhands, length(f), length(channels), length(data) - length(subjects_to_remove));

for i = 1 : (length(data) - length(subjects_to_remove))
    subj_name = data(i);

    % Check that the name is not associated with a non representative subject
    if ~contains(subj_name, 'aj7', 'IgnoreCase', true) && ~contains(subj_name, 'aj9', 'IgnoreCase', true) && ~contains(subj_name, 'ai7', 'IgnoreCase', true)
        ERDfeet_tot(:, :, :, i)    = subjects.(subj_name).ERDavg_feet(1:minLen_ERDfeet,:,:);
        ERDhands_tot(:, :, :, i)   = subjects.(subj_name).ERDavg_hands(1:minLen_ERDhands,:,:);
    end
end

% Compute Grand Average for representative subjects
ERDfeet_GA  = mean(ERDfeet_tot, 4);
ERDhands_GA = mean(ERDhands_tot, 4);


hf = figure('units','normalized','outerposition',[0 0 1 1]); 
hf.Name = 'GA ERD on PSD';
hf.NumberTitle = 'off';
subplot(231)
imagesc(t, f, ERDhands_GA(:, :, chns(1))')
colormap hot
colorbar
title('GA ERD avg both hands - C3')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(232)
imagesc(t, f, ERDhands_GA(:, :, chns(2))')
colormap hot
colorbar
title('GA ERD avg both hands - Cz')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(233)
imagesc(t, f, ERDhands_GA(:, :, chns(3))')
colormap hot
colorbar
title('GA ERD avg both hands - C4')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(234)
imagesc(t, f, ERDfeet_GA(:, :, chns(1))')
colormap hot
colorbar
title('GA ERD avg both feet - C3')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(235)
imagesc(t, f, ERDfeet_GA(:, :, chns(2))')
colormap hot
colorbar
title('GA ERD avg both feet - Cz')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(236)
imagesc(t, f, ERDfeet_GA(:, :, chns(3))')
colormap hot
colorbar
title('GA ERD avg both feet - C4')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')
