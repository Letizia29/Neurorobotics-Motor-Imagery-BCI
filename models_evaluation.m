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

% Data path
addpath(genpath(fullfile(pwd, "Data/")))

%% Data loading and storing

% As for the previous script, all the subjects' data are loaded and saved 
% in the 'subjects' struct

m = "_micontinuous";
data = [strcat("aj1",m),strcat("aj3",m),strcat("aj4",m),strcat("aj7",m),strcat("aj9",m),strcat("ai6",m),strcat("ai7",m),strcat("ai8",m)];

subjects = struct(); 

for i = 1:length(data)

    disp(['Loading subject ', data{i}(1:3), ' data'])
    subj_name = data(i);

    % Load data
    subjects.(subj_name) = load(fullfile(pwd, strcat("Data/", subj_name, '/data.mat')));
end

clc
disp('Done')


%% Features extraction

% From each subject's concatenated PSD file, the windows corresponding to
% the 3 most discriminative features are extracted (calculated from the 
% Fisher Score maps on the offline data)

for i = 1:length(data)
    subj_name = data(i);

    disp(['Extracting features for subject ', data{i}(1:3)])

    PSD_selected = zeros(size(subjects.(subj_name).data_subjects.wnds_CktoCFk_online, 1), length(subjects.(subj_name).data_subjects.col_feat));

    for j = 1:3
        % Row and column indeces of selected features
        row_feat = subjects.(subj_name).data_subjects.row_feat(j);
        col_feat = subjects.(subj_name).data_subjects.col_feat(j);
        
        PSD_selected(:, j) = subjects.(subj_name).data_subjects.wnds_CktoCFk_online(:, row_feat, col_feat);
    end
    
    % Save the value in the struct
    subjects.(subj_name).data_subjects.PSD_selected = PSD_selected;

end

clc
disp('Done')

%% Model prediction on online runs (single sample)

disp('Model evaluation on online runs')

for i = 1:length(data)
    subj_name = data(i);
    
    % Load variables
    model = subjects.(subj_name).mdl;
    PSD_selected = subjects.(subj_name).data_subjects.PSD_selected;
    Ck = subjects.(subj_name).data_subjects.class;

    % Overall training accuracy
    [Gk, pp] = predict(model, PSD_selected);
    
    % Both hands training accuracy
    [Gk_hands, pp_hands] = predict(model, PSD_selected(Ck == 773,:));

    % Both feet training accuracy
    [Gk_feet, pp_feet]   = predict(model, PSD_selected(Ck == 771,:));
    
    % Bar plot of single sample accuracies

    % Single sample accuracies
    overall_accuracy    = mean(Gk == Ck) * 100;                     % Overall accuracy in %
    both_hands_accuracy = mean(Gk_hands == Ck(Ck == 773)) * 100;    % Both hands accuracy in %
    both_feet_accuracy  = mean(Gk_feet == Ck(Ck == 771)) * 100;     % Both feet accuracy in %
    
    % Save variables in the data structure
    subjects.(subj_name).data_subjects.Gk = Gk;
    subjects.(subj_name).data_subjects.pp = pp;
    subjects.(subj_name).data_subjects.Gk_hands = Gk_hands;
    subjects.(subj_name).data_subjects.Gk_feet = Gk_feet;
    subjects.(subj_name).data_subjects.pp_hands = pp_hands;
    subjects.(subj_name).data_subjects.pp_feet = pp_feet;
    subjects.(subj_name).data_subjects.overall_acc = overall_accuracy;
    subjects.(subj_name).data_subjects.hands_acc = both_hands_accuracy;
    subjects.(subj_name).data_subjects.feet_acc = both_feet_accuracy;

    % Print accuracies
    fprintf('Accuracies of the model for subject %s\n', data{i}(1:3));
    fprintf('Accuracy: %f%%\n', overall_accuracy);
    fprintf('Accuracy both feet: %f%%\n', both_feet_accuracy);
    fprintf('Accuracy both hands: %f%%\n\n', both_hands_accuracy);

    % Data for the bar graph
    accuracies = [overall_accuracy, both_hands_accuracy, both_feet_accuracy];
    x_labels = {'Overall', 'Both Hands', 'Both Feet'};
    
    % Accuracies visualization as a bar graph
    hf = figure;
    hf.Name = ['Subject ', data{i}(1:3), ' accuracies on test'];
    hf.NumberTitle = 'off';
    bar(accuracies)
    set(gca, 'xticklabel', x_labels)
    ylabel('Accuracy [%]')
    title(['Single sample accuracy on test set - subject ', data{i}(1:3)])
    ylim([0, 100])
    grid on

end

%% Evidence accumulation framework

% Exponential accumulation framework
alpha = 0.95;                       % Smoothing parameter

for i = 1:length(data)
    subj_name = data(i);

    % Load posterior probabilities obtained from the prediction
    pp = subjects.(subj_name).data_subjects.pp;
    
    % For each window in Cue and Continuous Feedback tasks, extract the number of the trial
    trials_windows = subjects.(subj_name).data_subjects.vectors_PSD_online.Tk(subjects.(subj_name).data_subjects.vectors_PSD_online.Ak > 0 | subjects.(subj_name).data_subjects.vectors_PSD_online.CFk > 0);
    
    nwindows = length(pp);

    % Initialize Decision structure
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
    subjects.(subj_name).data_subjects.D = D;

end


%% Plot trial accuracy

% Define decision thresholds
thr = [0.3 0.7];

% Initialize data structure to store the average time to deliver a command 
t_avg = zeros(1, length(data));

for i = 1:length(data)
    subj_name = data(i);

    % For each window in Cue and Continuous Feedback tasks, extract the number of the trial
    trials_windows = subjects.(subj_name).data_subjects.vectors_PSD_online.Tk(subjects.(subj_name).data_subjects.vectors_PSD_online.Ak > 0 | subjects.(subj_name).data_subjects.vectors_PSD_online.CFk > 0);
    pp = subjects.(subj_name).data_subjects.pp;
    D = subjects.(subj_name).data_subjects.D;

    % Initialize structures
    Gk_trial_all = zeros(max(trials_windows),1);                % Store trial prediction for each trial
    t_for_command = zeros(max(trials_windows),1);               % Store the time for a command for each trial

    % Compute classification for succesful trials
    for trial_number = 1 : max(trials_windows)
        
        D_trial = D(trials_windows == trial_number,:);
        
        % Structure to store the sample when the threshold is reached
        t = zeros(1, length(D_trial));

        for j = 1 : length(D_trial)
            if D_trial(j) <= thr(1)
                Gk_trial_all(trial_number) = 773;           % Trial classified as both hands
                t(j) = j;                                   % Save the sample index if the th is reached
            end
            if D_trial(j) >= thr(2)
                Gk_trial_all(trial_number) = 771;           % Trial classified as both hands
                t(j) = j;                                   % Save the sample index if the th is reached
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
trial_accuracy_no_rejection = zeros(1, length(data));           % One value for each subject
trial_accuracy_no_rejection_feet  = zeros(1, length(data));
trial_accuracy_no_rejection_hands = zeros(1, length(data));

for i = 1:length(data)
    subj_name = data(i);

    % Store values in temporary variables
    trials_windows = subjects.(subj_name).data_subjects.vectors_PSD_online.Tk(subjects.(subj_name).data_subjects.vectors_PSD_online.Ak > 0 | subjects.(subj_name).data_subjects.vectors_PSD_online.CFk > 0);
    pp = subjects.(subj_name).data_subjects.pp;
    D = subjects.(subj_name).data_subjects.D;
    Gk_trial_all = subjects.(subj_name).Gk_trial_all;

    % Calculate accuracy
    trial_accuracy_no_rejection(i) = mean(Gk_trial_all == subjects.(subj_name).data_subjects.vectors_PSD_online.Ck) * 100;
    trial_accuracy_no_rejection_feet(i)  = mean(Gk_trial_all(subjects.(subj_name).data_subjects.vectors_PSD_online.Ck == 771) == subjects.(subj_name).data_subjects.vectors_PSD_online.Ck(subjects.(subj_name).data_subjects.vectors_PSD_online.Ck == 771)) * 100;
    trial_accuracy_no_rejection_hands(i) = mean(Gk_trial_all(subjects.(subj_name).data_subjects.vectors_PSD_online.Ck == 773) == subjects.(subj_name).data_subjects.vectors_PSD_online.Ck(subjects.(subj_name).data_subjects.vectors_PSD_online.Ck == 773)) * 100;

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
    trials_windows = subjects.(subj_name).data_subjects.vectors_PSD_online.Tk(subjects.(subj_name).data_subjects.vectors_PSD_online.Ak > 0 | subjects.(subj_name).data_subjects.vectors_PSD_online.CFk > 0);
    pp = subjects.(subj_name).data_subjects.pp;
    D = subjects.(subj_name).data_subjects.D;
    Gk_trial_all = subjects.(subj_name).Gk_trial_all;

    % Find the trials to reject
    idx_rejection = (Gk_trial_all ~= 0);

    % Calculate the accuracy on the remaining trials
    trial_accuracy_rejection(i) = mean(Gk_trial_all(idx_rejection) == subjects.(subj_name).data_subjects.vectors_PSD_online.Ck(idx_rejection)) * 100;
    trial_accuracy_rejection_feet(i)  = mean(Gk_trial_all(idx_rejection & subjects.(subj_name).data_subjects.vectors_PSD_online.Ck == 771) == subjects.(subj_name).data_subjects.vectors_PSD_online.Ck(idx_rejection & subjects.(subj_name).data_subjects.vectors_PSD_online.Ck == 771)) * 100;
    trial_accuracy_rejection_hands(i) = mean(Gk_trial_all(idx_rejection & subjects.(subj_name).data_subjects.vectors_PSD_online.Ck == 773) == subjects.(subj_name).data_subjects.vectors_PSD_online.Ck(idx_rejection & subjects.(subj_name).data_subjects.vectors_PSD_online.Ck == 773)) * 100;

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
    title(['Trial accuracy on test set with rejection - subject ', data{i}(1:3)])
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

%%

%trial_accuracy_rejection - trial_accuracy_no_rejection


% 
% % trial based accuracy
% % without rejection
% % accuracy = correct/total * 100
% trial_accuracy_no_rejection = mean(Gk_trial_all == vectors.Pk) * 100;
% 
% % remove rejected trials
% Gk_trial_with_rejection = Gk_trial_all(Gk_trial_all ~= 0);
% Pk_with_rejection = vectors.Pk(Gk_trial_all ~= 0);
% 
% trial_accuracy_rejection = mean(Gk_trial_with_rejection == Pk_with_rejection) * 100;
% 
% % Data for the bar graph
% trial_accuracies = [trial_accuracy_no_rejection, trial_accuracy_rejection];
% 
% % Labels for the x-axis
% x_labels = {'without rejection', 'with rejection'};
% 
% % Create the bar graph
% figure(4)
% bar(trial_accuracies)
% 
% % Customize the graph
% set(gca, 'xticklabel', x_labels)
% ylabel('Accuracy [%]')
% title('Trial accuracy on test set (with rest)')
% 
% grid on
% 
% 
% %% Trial accuracy with and without rejection of time-out trials (NO REST)
% Gk_trial_all_NR = Gk_trial_all(find(vectors.Pk ~= 783));
% Pk_NR = vectors.Pk(find(vectors.Pk ~= 783));
% Gk_trial_with_rejection_NR = Gk_trial_with_rejection(find(Pk_with_rejection ~= 783));
% Pk_with_rejection_NR = Pk_with_rejection(find(Pk_with_rejection ~= 783));
% 
% % trial based accuracy
% % accuracy = correct/total * 100
% 
% trial_accuracy_no_rejection_NR = mean(Gk_trial_all_NR == Pk_NR) * 100;
% trial_accuracy_rejection_NR = mean(Gk_trial_with_rejection_NR == Pk_with_rejection_NR) * 100;
% 
% % Data for the bar graph
% trial_accuracies = [trial_accuracy_no_rejection_NR, trial_accuracy_rejection_NR];
% 
% % Labels for the x-axis
% x_labels = {'without rejection', 'with rejection'};
% 
% % Create the bar graph
% figure(5)
% bar(trial_accuracies)
% 
% % Customize the graph
% set(gca, 'xticklabel', x_labels)
% ylabel('Accuracy [%]')
% title('Trial accuracy on test set (NO rest)')
% 
% grid on
