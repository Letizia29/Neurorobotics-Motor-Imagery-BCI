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

%% Data loading and storing

m = "_micontinuous";
data = [strcat("aj1",m),strcat("aj3",m),strcat("aj4",m),strcat("aj7",m),strcat("aj9",m),strcat("ai6",m),strcat("ai7",m),strcat("ai8",m)];

subjects = struct(); % where s and h data will be saved for each subject

for i = 1:length(data)

    disp(['Loading subject ', data{i}(1:3), ' data'])
    subj_name = data(i);

    % load data
    subjects.(subj_name) = load(fullfile(pwd, strcat("Data/", subj_name, '/data.mat')));
end

clc
disp('Done')


%% Features extraction

for i = 1:length(data)
    subj_name = data(i);

    PSD_selected = zeros(size(subjects.(subj_name).data_subjects.wnds_CktoCFk_online, 1), length(subjects.(subj_name).data_subjects.col_feat));

    for j = 1:3
        % row and column indeces of selected features
        row_ind = subjects.(subj_name).data_subjects.row_feat(j);
        col_ind = subjects.(subj_name).data_subjects.col_feat(j);
        
        PSD_selected(:, j) = subjects.(subj_name).data_subjects.wnds_CktoCFk_online(:, row_ind, col_ind);
    end
    
    subjects.(subj_name).data_subjects.PSD_selected = PSD_selected;

end

%% Model prediction on online runs

clc
disp('Model evaluation on online runs')

for i = 1:length(data)
    subj_name = data(i);
    
    % necessary variables
    model = subjects.(subj_name).mdl;
    PSD_selected = subjects.(subj_name).data_subjects.PSD_selected;
    Ck = subjects.(subj_name).data_subjects.class;

    % overall training accuracy
    [Gk, pp] = predict(model, PSD_selected);
    
    % both hands training accuracy
    [Gk_hands, pp_hands] = predict(model, PSD_selected(Ck == 773,:));
    % both feet training accuracy
    [Gk_feet, pp_feet] = predict(model, PSD_selected(Ck == 771,:));
    
    % bar plot of single sample accuracies
    % Sample accuracies
    overall_accuracy    = mean(Gk == Ck) * 100;                     % Overall accuracy in %
    both_hands_accuracy = mean(Gk_hands == Ck(Ck == 773)) * 100;    % Both hands accuracy in %
    both_feet_accuracy  = mean(Gk_feet == Ck(Ck == 771)) * 100;     % Both feet accuracy in %
    
    % save in the data structure
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
    fprintf('Accuracy: %f\n', overall_accuracy);
    fprintf('Accuracy both feet: %f\n', both_hands_accuracy);
    fprintf('Accuracy both hands: %f\n\n', both_feet_accuracy);

    % Data for the bar graph
    accuracies = [overall_accuracy, both_hands_accuracy, both_feet_accuracy];
    x_labels = {'Overall', 'Both Hands', 'Both Feet'};
    
    figure(i)
    bar(accuracies)
    set(gca, 'xticklabel', x_labels)
    ylabel('Accuracy [%]')
    title('Single sample accuracy on test set')
    ylim([0, 100])
    grid on

end

%% Control framework
% (trial based accuracy)

% Exponential accumulation framework
alpha = 0.95; % smoothing parameter [0 1]

for i = 1:length(data)
    subj_name = data(i);
    pp = subjects.(subj_name).data_subjects.pp;

    % COME SCUSA
    trials_windows = subjects.(subj_name).data_subjects.vectors_PSD_online.Tk(subjects.(subj_name).data_subjects.vectors_PSD_online.Ak > 0 | subjects.(subj_name).data_subjects.vectors_PSD_online.CFk > 0);
    
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
    subjects.(subj_name).data_subjects.D = D;

end


% DA ADATTARE

%% Plot trial accuracy

thr = [0.2 0.8];

% posterior proabablity pp of trial 55
% trial_number = randi(100, 1);
sel_trial = 30;

for i = 1:length(data)
    subj_name = data(i);

    trials_windows = subjects.(subj_name).data_subjects.vectors_PSD_online.Tk(subjects.(subj_name).data_subjects.vectors_PSD_online.Ak > 0 | subjects.(subj_name).data_subjects.vectors_PSD_online.CFk > 0);
    pp = subjects.(subj_name).data_subjects.pp;
    D = subjects.(subj_name).data_subjects.D;

    pp_trial = pp(trials_windows == sel_trial,:);
    D_trial = D(trials_windows == sel_trial,:); % integrated probability
    samples = 1:length(pp_trial(:,2));

    figure
    hold on
    scatter(samples, pp_trial(:,2), 'k')
    plot(samples, D_trial(:,2), 'k', 'LineWidth', 2)
    title(['Trial ', num2str(sel_trial), ' - Class both hands'])
    xlabel('samples')
    ylabel('probability/control')
    legend('posterior probability', 'integrated probability', 'Location', 'best')
    yline(0.5, '--')
    yline(thr(1), 'k')
    yline(thr(2), 'k')
    ylim([0 1])
    xlim([1 samples(end)])

    Gk_trial_all = zeros(max(trials_windows),1);

    % compute classification for succesful trials (no timeout)
    for trial_number = 1 : max(trials_windows)
        D_trial = D(trials_windows == trial_number,:);
        for j = 1 : length(D_trial)
            if D_trial(j) <= thr(1)
                Gk_trial_all(trial_number) = 773; % trial classified as both hands
                % break ?????????????????????
            end
            if D_trial(j) >= thr(2)
                Gk_trial_all(trial_number) = 771; % trial classified as both hands
                % break ?????????????????????
                % Metterei break (se ho capito bene cosa fa il codice)
                % in modo tale da aggiornare il vettore Gk_trial_all e
                % uscire subito, senza poterlo più aggiornare dopo, appena
                % si supera la soglia

                % La mia implementazione (non penso sia giusta ma avevo pensato di fare sta cosa)
                % era quella di contare il numero di sample che superano la
                % soglia e prendere la classe relativa al contatore con
                % valore maggiore. in questo modo, anche se per sbaglio si
                % supera la soglia sbagliata, se si riesce a correggere il
                % numero di campioni contati "aggiusta" la predizione

                % Però a te veniva simile al prof, quindi direi di lasciare
                % la tua
            end
        end
    end

    subjects.(subj_name).Gk_trial_all = Gk_trial_all;
end



%% Trial accuracy with and without rejection of time-out trials

trial_accuracy_no_rejection = zeros(1, length(data));

for i = 1:length(data)
    subj_name = data(i);

    trials_windows = subjects.(subj_name).data_subjects.vectors_PSD_online.Tk(subjects.(subj_name).data_subjects.vectors_PSD_online.Ak > 0 | subjects.(subj_name).data_subjects.vectors_PSD_online.CFk > 0);
    pp = subjects.(subj_name).data_subjects.pp;
    D = subjects.(subj_name).data_subjects.D;

    

    Gk_trial_all = subjects.(subj_name).Gk_trial_all;

    trial_accuracy_no_rejection(i) = mean(Gk_trial_all == subjects.(subj_name).data_subjects.vectors_PSD_online.Ck) * 100;

end
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
