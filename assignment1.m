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

addpath(genpath(fullfile(pwd, "Data/")))

% Loading Laplacian Mask
load('laplacian16.mat');

%% Data loading and PSD computation
% for each subject, loading of each run, application of Laplacian filter
% and PSD computation, with consequent events re-computing

m = "_micontinuous";
data = [strcat("aj1",m),strcat("aj3",m),strcat("aj4",m),strcat("aj7",m),strcat("aj9",m),strcat("ai6",m),strcat("ai7",m),strcat("ai8",m)];

subjects = struct(); % where s and h data will be saved for each subject

disp('Loading subjects data')


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


                [b, a] = butter(3, [3/(512/2), 49/(512/2)]);
                s = filtfilt(b, a, s);

                % Remove artifacts
                [s, h] = remArtifacts(s, h);

                % t = 0 : 1/512 : (size(s, 1) - 1)/512;
                % figure()
                % plot(t, s)
                % title(['subj ', num2str(i), ' run ', num2str(j-2)])

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

sample_rate = 512; % [Hz]


%% Concatenate signals and events' types, positions, durations for each subject

clc

for i = 1:length(data)

    disp(['Subject ', data{i}(1:3)])

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

clc
disp('Done')

%% Data processing

% FILTERING
% Butterworth filters
n_mu = 4; % filter order
n_beta = 4;

% mu band
W1 = 8; % Hz
W2 = 12; % Hz
Wn_mu = 2*[W1 W2]/sample_rate;

% beta band
W1 = 16; % Hz
W2 = 24; % Hz
Wn_beta = 2*[W1 W2]/sample_rate;

% filter coefficients
[b_mu, a_mu] = butter(n_mu, Wn_mu);
[b_beta, a_beta] = butter(n_beta, Wn_beta);

starting_time_cue = zeros(length(data), 1);
starting_time_cf  = zeros(length(data), 1);

for i = 1:length(data)

    disp(['Processing subject ', data{i}(1:3), ' data'])

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
    B = ones(1, LengthWin*sample_rate)/LengthWin/sample_rate;
    % filter the signal
    sfilt_sq_ma_mu = filter(B, A, sfilt_sq_mu);
    sfilt_sq_ma_beta = filter(B, A, sfilt_sq_beta);
    
    % LOG BAND POWER
    % Logarithm transform
    logBP_mu = log(sfilt_sq_ma_mu);
    logBP_beta = log(sfilt_sq_ma_beta);

    % save results
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

    subjects.(subj_name).ERD_logBP_mu = 100 * (Activity_mu - Reference_mu) ./ Reference_mu;
    subjects.(subj_name).ERD_logBP_beta = 100 * (Activity_beta - Reference_beta) ./ Reference_beta;
    

end 

clc
disp('Done')


%% Plot ERD/ERS on logBP

% significant channels (C3, Cz, C4)
chns = [7, 9, 11];

ax1 = figure(1); 
ax1.Name = 'Average ERD logBP \mu - channel C3'; % 7
ax1.Position = [50,100,1600,600];
hold off

ax2 = figure(2); 
ax2.Name = 'Average ERD logBP \mu - channel Cz'; % 9
ax2.Position = [50,100,1600,600];
hold off

ax3 = figure(3); 
ax3.Name = 'Average ERD logBP \mu - channel C4'; % 11
ax3.Position = [50,100,1600,600];
hold off

ax4 = figure(4); 
ax4.Name = 'Average ERD logBP \beta - channel C3'; % 7
ax4.Position = [50,100,1600,600];
hold off

ax5 = figure(5); 
ax5.Name = 'Average ERD logBP \beta - channel Cz'; % 9
ax5.Position = [50,100,1600,600];
hold off

ax6 = figure(6); 
ax6.Name = 'Average ERD logBP \beta - channel C4'; % 11
ax6.Position = [50,100,1600,600];
hold off

for i = 1:length(data)
    subj_name = data(i);

    ERDmu_feet    = subjects.(subj_name).ERD_logBP_mu(:, :, subjects.(subj_name).vectors.Ck == 771);
    ERDmu_hands   = subjects.(subj_name).ERD_logBP_mu(:, :, subjects.(subj_name).vectors.Ck == 773);
    ERDbeta_feet  = subjects.(subj_name).ERD_logBP_beta(:, :, subjects.(subj_name).vectors.Ck == 771);
    ERDbeta_hands = subjects.(subj_name).ERD_logBP_beta(:, :, subjects.(subj_name).vectors.Ck == 773);

    % average ERD/ERS
    ERD_logBP_mu_avg_feet  = mean(ERDmu_feet, 3);
    ERD_logBP_mu_avg_hands = mean(ERDmu_hands, 3);

    ERD_logBP_beta_avg_feet  = mean(ERDbeta_feet, 3);
    ERD_logBP_beta_avg_hands = mean(ERDbeta_hands, 3);

    % standard deviation
    ERD_logBP_mu_SE_feet  = std(ERDmu_feet, 0, 3)./sqrt(length(ERDmu_feet(1,1,:)));
    ERD_logBP_mu_SE_hands = std(ERDmu_hands, 0, 3)./sqrt(length(ERDmu_hands(1,1,:)));

    ERD_logBP_beta_SE_feet  = std(ERDbeta_feet, 0, 3)./sqrt(length(ERDbeta_feet(1,1,:)));
    ERD_logBP_beta_SE_hands = std(ERDbeta_hands, 0, 3)./sqrt(length(ERDbeta_hands(1,1,:)));


    % Visualization
    

    % time vector

    % ????????
    T = 0.0627; % wshift

    T = 1/sample_rate;
    t = 0:T:(length(ERD_logBP_mu_avg_feet)-1)*T;

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

%% COMMENTS ON THE PLOTS

% MU BAND

%   - Channel C3:
% There is a clear desynchronization of logBP during hands task for
% subjects 1, 3 and 8
% The ERD is also present for subjects 2, 4, 5 and 6. However, the ERD is
% not very clear and it is followed by high ERS

%   - Channel Cz:
% 

%   - Channel C4:
%

% BETA BAND

%   - Channel C3:
%

%   - Channel Cz:
%

%   - Channel C4:
%


%% GRAND AVERAGE 
% for ERD on logarithmic band power (except sub 7)

% adjust different lengths
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

% save data from all subjects
ERDmu_feet_tot    = zeros(minLen_ERDmu_feet, 16, length(data)-1);
ERDmu_hands_tot   = zeros(minLen_ERDmu_hands, 16, length(data)-1);
ERDbeta_feet_tot  = zeros(minLen_ERDbeta_feet, 16, length(data)-1);
ERDbeta_hands_tot = zeros(minLen_ERDbeta_hands, 16, length(data)-1);

remove_sub7 = 0;

if remove_sub7 == 1

    for i = 1:length(data)
        subj_name = data(i);
        if i < 7 % remove subject 7
            ERDmu_feet_tot(:, :, i)    = subjects.(subj_name).ERD_logBP_mu_avg_feet(1:minLen_ERDmu_feet,:);
            ERDmu_hands_tot(:, :, i)   = subjects.(subj_name).ERD_logBP_mu_avg_hands(1:minLen_ERDmu_hands,:);
            ERDbeta_feet_tot(:, :, i)  = subjects.(subj_name).ERD_logBP_beta_avg_feet(1:minLen_ERDbeta_feet,:);
            ERDbeta_hands_tot(:, :, i) = subjects.(subj_name).ERD_logBP_beta_avg_hands(1:minLen_ERDbeta_hands,:);
        end
        if i == 8
            ERDmu_feet_tot(:, :, i-1)    = subjects.(subj_name).ERD_logBP_mu_avg_feet(1:minLen_ERDmu_feet,:);
            ERDmu_hands_tot(:, :, i-1)   = subjects.(subj_name).ERD_logBP_mu_avg_hands(1:minLen_ERDmu_hands,:);
            ERDbeta_feet_tot(:, :, i-1)  = subjects.(subj_name).ERD_logBP_beta_avg_feet(1:minLen_ERDbeta_feet,:);
            ERDbeta_hands_tot(:, :, i-1) = subjects.(subj_name).ERD_logBP_beta_avg_hands(1:minLen_ERDbeta_hands,:);
        end
    end
else
    for i = 1:length(data)
        subj_name = data(i);
        ERDmu_feet_tot(:, :, i)    = subjects.(subj_name).ERD_logBP_mu_avg_feet(1:minLen_ERDmu_feet,:);
        ERDmu_hands_tot(:, :, i)   = subjects.(subj_name).ERD_logBP_mu_avg_hands(1:minLen_ERDmu_hands,:);
        ERDbeta_feet_tot(:, :, i)  = subjects.(subj_name).ERD_logBP_beta_avg_feet(1:minLen_ERDbeta_feet,:);
        ERDbeta_hands_tot(:, :, i) = subjects.(subj_name).ERD_logBP_beta_avg_hands(1:minLen_ERDbeta_hands,:);

    end
end

% compute Grand Average
ERDmu_feet_GA    = mean(ERDmu_feet_tot, 3);
ERDmu_hands_GA   = mean(ERDmu_hands_tot, 3);
ERDbeta_feet_GA  = mean(ERDbeta_feet_tot, 3);
ERDbeta_hands_GA = mean(ERDbeta_hands_tot, 3);

ERDmu_feet_SE    = std(ERDmu_feet_tot, 0, 3)./sqrt(length(ERDmu_feet_tot(1, 1, :)));
ERDmu_hands_SE   = std(ERDmu_hands_tot, 0, 3)./sqrt(length(ERDmu_hands_tot(1, 1, :)));
ERDbeta_feet_SE  = std(ERDbeta_feet_tot, 0, 3)./sqrt(length(ERDbeta_feet_tot(1, 1, :)));
ERDbeta_hands_SE = std(ERDbeta_hands_tot, 0, 3)./sqrt(length(ERDbeta_hands_tot(1, 1, :)));

% GA temporal plots
figure(7)
% time axis 
T = 1/sample_rate;
t = 0:T:(length(ERDmu_feet_GA(:,1,1))-1)*T;

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
% ylim([-60, 170])
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
% ylim([-60, 170])
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
% ylim([-60, 170])
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
% ylim([-60, 170])
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
% ylim([-60, 170])
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
% ylim([-60, 170])
title('Grand Average ERD logBP \beta C4')
legend('both feet', 'both hands')
hold off


% Topographic plots
load('chanlocs16.mat');
len = min(h.DUR(h.TYP == 786)); % presa dall'ultima h salvata giusto per fare un prova

% mu band
ERD_Ref_773 = mean(ERDmu_hands_GA(1:len, :), 1);
ERD_Act_773 = mean(ERDmu_hands_GA(len+1:end, :), 1);
ERD_Ref_771 = mean(ERDmu_feet_GA(1:len, :), 1);
ERD_Act_771 = mean(ERDmu_feet_GA(len+1:end, :), 1);

figure(8);
subplot(221)
topoplot(squeeze(ERD_Ref_773), chanlocs16);
title('Reference - \mu band - both hands')
colorbar
clim([-20, 50])
subplot(222)
topoplot(squeeze(ERD_Act_773), chanlocs16);
title('Activity - \mu band - both hands')
colorbar
clim([-20, 50])
subplot(223)
topoplot(squeeze(ERD_Ref_771), chanlocs16);
title('Reference - \mu band - both feet')
colorbar
clim([-20, 50])
subplot(224)
topoplot(squeeze(ERD_Act_771), chanlocs16);
title('Activity - \mu band - both feet')
colorbar
clim([-20, 50])

% beta band
ERD_Ref_773 = mean(ERDbeta_hands_GA(1:len, :), 1);
ERD_Act_773 = mean(ERDbeta_hands_GA(len+1:end, :), 1);
ERD_Ref_771 = mean(ERDbeta_feet_GA(1:len, :), 1);
ERD_Act_771 = mean(ERDbeta_feet_GA(len+1:end, :), 1);

figure(9);
subplot(221)
topoplot(squeeze(ERD_Ref_773), chanlocs16);
title('Reference - \beta band - both hands')
colorbar
clim([-20, 50])
subplot(222)
topoplot(squeeze(ERD_Act_773), chanlocs16);
title('Activity - \beta band - both hands')
colorbar
clim([-20, 50])
subplot(223)
topoplot(squeeze(ERD_Ref_771), chanlocs16);
title('Reference - \beta band - both feet')
colorbar
clim([-20, 50])
subplot(224)
topoplot(squeeze(ERD_Act_771), chanlocs16);
title('Activity - \beta band - both feet')
colorbar
clim([-20, 50])


%% ERD/ERS Spectrogram

%% Concatenate the files

for i = 1:length(data)

    disp(['Subject ', data{i}(1:3)])

    subj_name = data(i);
    runs_names = fieldnames(subjects.(subj_name));
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

    subjects.(subj_name).h_PSD.POS = POS;
    subjects.(subj_name).h_PSD.DUR = DUR;
    subjects.(subj_name).h_PSD.TYP = subjects.(subj_name).h.TYP;

    subjects.(subj_name).vectors_PSD = labelVecs(subjects.(subj_name).PSD_c, subjects.(subj_name).h_PSD);

    % Concatenation of the PSD for the online files

    subjects.(subj_name).PSD_c_online = [];
    POS = [];
    DUR = [];

    for j = 1:length(runs_names)

        if contains(runs_names{j}, 'online', 'IgnoreCase', true)

            DUR = [DUR; subjects.(subj_name).(runs_names{j}).h_PSD.DUR];
            POS = [POS; subjects.(subj_name).(runs_names{j}).h_PSD.POS + length(subjects.(subj_name).PSD_c_online)];

            subjects.(subj_name).PSD_c_online = [subjects.(subj_name).PSD_c_online; subjects.(subj_name).(runs_names{j}).PSD];        
        end

    end

    subjects.(subj_name).h_PSD_online.POS = POS;
    subjects.(subj_name).h_PSD_online.DUR = DUR;
    subjects.(subj_name).h_PSD_online.TYP = subjects.(subj_name).h.TYP;

    subjects.(subj_name).vectors_PSD_online = labelVecs(subjects.(subj_name).PSD_c_online, subjects.(subj_name).h_PSD_online);
    
end

%% Activity and Reference computation

ax1 = figure(10); 
ax1.Name = 'Average ERD for both feet task - channel C3'; % 7
ax1.Position = [50,100,1600,600];
hold off

ax2 = figure(11); 
ax2.Name = 'Average ERD for both feet task - channel Cz'; % 9
ax2.Position = [50,100,1600,600];
hold off

ax3 = figure(12); 
ax3.Name = 'Average ERD for both feet task - channel C4'; % 11
ax3.Position = [50,100,1600,600];
hold off

ax4 = figure(13); 
ax4.Name = 'Average ERD for both hands task - channel C3'; % 7
ax4.Position = [50,100,1600,600];
hold off

ax5 = figure(14); 
ax5.Name = 'Average ERD for both hands task - channel Cz'; % 9
ax5.Position = [50,100,1600,600];
hold off

ax6 = figure(15); 
ax6.Name = 'Average ERD for both hands task - channel C4'; % 11
ax6.Position = [50,100,1600,600];
hold off

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

    ERD = 100 * (Activity- Reference)./ Reference;
    % ERD = log(Activity ./ Reference);
    % ERD = 10 * log10(Activity ./ Reference);
    
    subjects.(subj_name).ERD_PSD = ERD;

    % ho aggiunto la definizione di Ck dentro labelVecs, perché l'ho usato anche sopra,
    % si può prendere da là direttamente
    % Ck = h_PSD.TYP(h_PSD.TYP == 771 | h_PSD.TYP == 773);

    ERDavg_feet  = mean(ERD(:, :, :, subjects.(subj_name).vectors_PSD.Ck == 771), 4);
    ERDavg_hands = mean(ERD(:, :, :, subjects.(subj_name).vectors_PSD.Ck == 773), 4);

    subjects.(subj_name).ERDavg_feet = ERDavg_feet;
    subjects.(subj_name).ERDavg_hands = ERDavg_hands;


    % Visualization
    
    % frequency vector
    f = subjects.(subj_name).(runs_names{1}).f;
    % time vector
    T = 0.0627; % wshift
    t = 0:T:length(ERDavg_hands)*T;

    % FARE QUALCOSA SU QUESTE OSCENITA'

    set(0,'CurrentFigure',ax1)
    subplot(2, 4, mod(i-1, 8)+1)
    imagesc(t, f, ERDavg_feet(:, :, 7)')
    colormap hot
    colorbar
    % clim([-1.1, 1.7])
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
    % clim([-1.1, 1.7])
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
    % clim([-1.1, 1.7])
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
    % clim([-1.1, 1.7])
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
    % clim([-1.1, 1.7])
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
    % clim([-1.1, 1.7])
    name = char(subj_name);
    title(strcat('ERD avg both hands - C4 - ',name(1:3)))
    set(gca,'YDir','normal')
    drawnow
    hold off

end

%% GRAND AVERAGE
% for ERS/ERS on PSD

% adjust different lengths
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

% save data from all subjects
ERDfeet_tot    = zeros(minLen_ERDfeet, 23, 16, length(data)-1);
ERDhands_tot   = zeros(minLen_ERDhands, 23, 16, length(data)-1);

remove_sub7 = 0;

if remove_sub7 == 1

    for i = 1:length(data)
        subj_name = data(i);
        if i < 7 % remove subject 7
            ERDfeet_tot(:, :, :, i)    = subjects.(subj_name).ERDavg_feet(1:minLen_ERDfeet,:,:);
            ERDhands_tot(:, :, :, i)   = subjects.(subj_name).ERDavg_hands(1:minLen_ERDhands,:,:);
        end
        if i == 8
            ERDfeet_tot(:, :, :, i-1)    = subjects.(subj_name).ERDavg_feet(1:minLen_ERDfeet,:,:);
            ERDhands_tot(:, :, :, i-1)   = subjects.(subj_name).ERDavg_hands(1:minLen_ERDhands,:,:);
        end
    end
else
    for i = 1:length(data)
        subj_name = data(i);
        ERDfeet_tot(:, :, :, i)    = subjects.(subj_name).ERDavg_feet(1:minLen_ERDfeet,:,:);
        ERDhands_tot(:, :, :, i)   = subjects.(subj_name).ERDavg_hands(1:minLen_ERDhands,:,:);
    end
end

% compute Grand Average
ERDfeet_GA    = mean(ERDfeet_tot, 4);
ERDhands_GA   = mean(ERDhands_tot, 4);


figure(30)

subplot(231)
imagesc(t, f, ERDhands_GA(:, :, chns(1))')
colormap hot
colorbar
% clim([-1.1, 1.7])
name = char(subj_name);
title(strcat('ERD avg both hands - C3 - ',name(1:3)))
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(232)
imagesc(t, f, ERDhands_GA(:, :, chns(2))')
colormap hot
colorbar
% clim([-1.1, 1.7])
name = char(subj_name);
title(strcat('ERD avg both hands - Cz - ',name(1:3)))
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(233)
imagesc(t, f, ERDhands_GA(:, :, chns(3))')
colormap hot
colorbar
% clim([-1.1, 1.7])
name = char(subj_name);
title(strcat('ERD avg both hands - C4 - ',name(1:3)))
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(234)
imagesc(t, f, ERDfeet_GA(:, :, chns(1))')
colormap hot
colorbar
% clim([-1.1, 1.7])
name = char(subj_name);
title(strcat('ERD avg both feet - C3 - ',name(1:3)))
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(235)
imagesc(t, f, ERDfeet_GA(:, :, chns(2))')
colormap hot
colorbar
% clim([-1.1, 1.7])
name = char(subj_name);
title(strcat('ERD avg both feet - Cz - ',name(1:3)))
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')

subplot(236)
imagesc(t, f, ERDfeet_GA(:, :, chns(3))')
colormap hot
colorbar
% clim([-1.1, 1.7])
name = char(subj_name);
title(strcat('ERD avg both feet - C4 - ',name(1:3)))
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca,'YDir','normal')


%% Feature selection

% NB: SI PUO' METTERE NEL CICLO PRIMA MA PER CHIAREZZA LO FACCIO FUORI. POI
% IN CASO BASTA UNIRLI.

channels = {"Fz", "FC3", "FC1", "FCz", "FC2", "FC4", "C3", "C1", "Cz", "C2", "C4", "CP3", "CP1", "CPz", "CP2", "CP4"};

ax0 = figure(16);

for i = 1:length(data)
    subj_name = data(i);

    % Define features
    wnds_CktoCFk = log(subjects.(subj_name).PSD_c(subjects.(subj_name).vectors_PSD.Ak > 0 | subjects.(subj_name).vectors_PSD.CFk > 0, :, :));
    features = reshape(wnds_CktoCFk, [size(wnds_CktoCFk, 1), size(wnds_CktoCFk, 2) * size(wnds_CktoCFk, 3)]);

    idx = subjects.(subj_name).vectors_PSD.Ak + subjects.(subj_name).vectors_PSD.CFk;       % Vector containing 771, 773 or 781
    idx = idx(idx > 0);
    class = (idx < 781) .* idx;     % Vector containing the class for each window

    for j = 2 : length(idx)
        if class(j-1) > 0 && class(j) == 0      % If the current class is 0, then
            class(j) = class(j-1);              % the value is updated according
        end                                     % to the previous state.
    end

    % Fisher score computation
    FS = abs(mean(features(class == 771, :), 1) - mean(features(class == 773, :), 1))./sqrt(std(features(class == 771, :), 1).^2 + std(features(class == 773, :), 1).^2);

    set(0,'CurrentFigure',ax0)
    subplot(2, 4, mod(i-1, 8)+1, 'Parent', ax0)
    title(strcat('Subject ', num2str(i)))
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

    idx = subjects.(subj_name).vectors_PSD_online.Ak + subjects.(subj_name).vectors_PSD_online.CFk;       % Vector containing 771, 773 or 781
    idx = idx(idx > 0);
    class = (idx ~= 781) .* idx;     % Vector containing the class for each window

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


% Identify and extract the most relevant features for each subject and on
% subject average


%% Models

for i = 1:length(data)
    subj_name = data(i);

    % Define the new matrix containing the data corresponding to the most
    % discriminative features

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

    data_subjects.PSD_c_online = subjects.(subj_name).PSD_c_online;
    data_subjects.h_PSD_online = subjects.(subj_name).h_PSD_online;
    data_subjects.vectors_PSD_online = subjects.(subj_name).vectors_PSD_online;
    data_subjects.wnds_CktoCFk_online = subjects.(subj_name).wnds_CktoCFk_online;
    data_subjects.features_online = subjects.(subj_name).features_online;
    data_subjects.row_feat = subjects.(subj_name).row_feat;
    data_subjects.col_feat = subjects.(subj_name).col_feat;
    data_subjects.class = subjects.(subj_name).class_online;

    name = strcat(pwd, '\Data\', subj_name, '\data');
    save(name, 'mdl', 'data_subjects')
    
end

%% Model evaluation on training data

for i = 1:length(data)
    subj_name = data(i);

    mdl = subjects.(subj_name).mdl;

    % Model prediction
    [Gk, pp] = predict(mdl, subjects.(subj_name).train_set);

    y = subjects.(subj_name).class;

    % Calculate the total accuracy (correct predictions over total predictions)
    tot_accuracy = mean(Gk == y) * 100;
    
    % Calculate class accuracy separately
    feet_accuracy  = mean(Gk(y == 771) == y(y == 771)) * 100;
    hands_accuracy = mean(Gk(y == 773) == y(y == 773)) * 100;
    
    % Print accuracies
    fprintf('Accuracies of the model for subject %s\n', data{i}(1:3));
    fprintf('Accuracy: %f\n', tot_accuracy);
    fprintf('Accuracy both feet: %f\n', feet_accuracy);
    fprintf('Accuracy both hands: %f\n\n', hands_accuracy);

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


%% Accumulation framework for online data

% (trial based accuracy)

% Exponential accumulation framework
alpha = 0.95; % smoothing parameter [0 1]

for i = 1:length(data)
    subj_name = data(i);
    pp = subjects.(subj_name).pp;

    % COME SCUSA
    trials_windows = subjects.(subj_name).vectors_PSD.Tk(subjects.(subj_name).vectors_PSD.Ak > 0 | subjects.(subj_name).vectors_PSD.CFk > 0);
    
    nwindows = length(pp);
    D = 0.5 * ones(nwindows, 2);
    for wId = 2:nwindows
        % Is the first sample of a new trial?
        if trials_windows(wId) ~= trials_windows(wId-1)
            % |- YES: Reset the current evidence (i.e., D(wId) to [0.5 0.5])
            D(wId,:) = [0.5 0.5];
        else 
            % |- NO:  Keep integrating the value
            D(wId,:) = D(wId - 1,:) * alpha + pp(wId,:) * (1 - alpha);
        end
    end
    
    % save results
    subjects.(subj_name).D = D;

end


% DA ADATTARE

%% Plot trial accuracy

thr = [0.2 0.8];

% posterior proabablity pp of trial 55
% trial_number = randi(100, 1);
sel_trial = 30;

t_avg = zeros(1, length(data));

for i = 1:length(data)
    subj_name = data(i);

    % salvare il tempo medio per raggiungere la threshold

    trials_windows = subjects.(subj_name).vectors_PSD.Tk(subjects.(subj_name).vectors_PSD.Ak > 0 | subjects.(subj_name).vectors_PSD.CFk > 0);
    pp = subjects.(subj_name).pp;
    D = subjects.(subj_name).D;

    pp_trial = pp(trials_windows == sel_trial,:);
    D_trial = D(trials_windows == sel_trial,:); % integrated probability
    samples = 1:length(pp_trial(:,2));

    % figure
    % hold on
    % scatter(samples, pp_trial(:,2), 'k')
    % plot(samples, D_trial(:,2), 'k', 'LineWidth', 2)
    % title(['Trial ', num2str(sel_trial), ' - Class both hands'])
    % xlabel('samples')
    % ylabel('probability/control')
    % legend('posterior probability', 'integrated probability', 'Location', 'best')
    % yline(0.5, '--')
    % yline(thr(1), 'k')
    % yline(thr(2), 'k')
    % ylim([0 1])
    % xlim([1 samples(end)])

    Gk_trial_all = zeros(max(trials_windows),1);

    t_for_command = zeros(max(trials_windows),1);

    % compute classification for succesful trials (no timeout)
    for trial_number = 1 : max(trials_windows)
        
        
        D_trial = D(trials_windows == trial_number,:);

        t = zeros(1, length(D_trial));

        for j = 1 : length(D_trial)
            if D_trial(j) <= thr(1)
                Gk_trial_all(trial_number) = 773; % trial classified as both hands
                t(j) = j;
                % breaks
            end
            if D_trial(j) >= thr(2)
                Gk_trial_all(trial_number) = 771; % trial classified as both hands
                t(j) = j;
                % break
                % contatore??????

            end
        end

        if sum(t) > 0
            t_for_command(trial_number) = find(t > 0, 1, "first");
        end

    end
    t_avg(i) = mean(t_for_command(t_for_command > 0));

    subjects.(subj_name).Gk_trial_all = Gk_trial_all;
end

% Dobbiamo trasformare in secondi
t_avg = t_avg * 0.0625;


%% Trial accuracy with and without rejection of time-out trials

trial_accuracy_no_rejection = zeros(1, length(data));
trial_accuracy_no_rejection_feet  = zeros(1, length(data));
trial_accuracy_no_rejection_hands = zeros(1, length(data));

for i = 1:length(data)
    subj_name = data(i);

    trials_windows = subjects.(subj_name).vectors_PSD.Tk(subjects.(subj_name).vectors_PSD.Ak > 0 | subjects.(subj_name).vectors_PSD.CFk > 0);
    pp = subjects.(subj_name).pp;
    D = subjects.(subj_name).D;

    

    Gk_trial_all = subjects.(subj_name).Gk_trial_all;

    trial_accuracy_no_rejection(i) = mean(Gk_trial_all == subjects.(subj_name).vectors_PSD.Ck) * 100;
    trial_accuracy_no_rejection_feet(i)  = mean(Gk_trial_all(subjects.(subj_name).vectors_PSD.Ck == 771) == subjects.(subj_name).vectors_PSD.Ck(subjects.(subj_name).vectors_PSD.Ck == 771)) * 100;
    trial_accuracy_no_rejection_hands(i) = mean(Gk_trial_all(subjects.(subj_name).vectors_PSD.Ck == 773) == subjects.(subj_name).vectors_PSD.Ck(subjects.(subj_name).vectors_PSD.Ck == 773)) * 100;

    accuracies = [trial_accuracy_no_rejection(i), trial_accuracy_no_rejection_feet(i), trial_accuracy_no_rejection_hands(i)];
    x_labels = {'Overall', 'Both Hands', 'Both Feet'};
    
    figure()
    bar(accuracies)
    set(gca, 'xticklabel', x_labels)
    ylabel('Accuracy [%]')
    title('Trial accuracy on test set')
    ylim([0, 100])
    grid on
end
