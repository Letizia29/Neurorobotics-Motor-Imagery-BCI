function [s, h] = remArtifacts(s, h)

% This function detects artifacts in raw data and removes the contaminated
% trials

diff_s = [zeros(1, 16); diff(s)];                   % Used for artifact detection

% Artifact condition
artifact_condition = (abs(s) > 150 & abs(diff_s) > 125);      % Valori da cambiare !!!!
[row, ~] = find(artifact_condition);

% Initialization
trial_flag = zeros(length(h.EVENT.POS)/4);          % Boolean information for each trial
h_bool = zeros(size(h.EVENT.POS));                  % Where the artifact is present in h.EVENT
s_bool = zeros(size(s, 1), 1);                      % Where the artifact is present s
trial_DUR = zeros(length(h.EVENT.POS)/4);           % Duration of each trial (to recompute POS at the end)

% Create the final values of the positions
h_POS = h.EVENT.POS;

% Find to which trial they correspond
for trialID = 1 : (length(h.EVENT.POS)/4)

    start_POS = h.EVENT.POS(1 + 4*(trialID-1));
    end_POS   = h.EVENT.POS(4*(trialID)) + h.EVENT.DUR(4*(trialID))-1;
    trial_DUR(trialID) = end_POS - start_POS;

    POS = start_POS : end_POS;                      % Vector of the positions
    if sum(ismember(row, POS)) > 0                  % Sample containing the artifact in POS
        trial_flag(trialID) = 1;                    % Mark the trial

        % Recompute the position vector
        h_POS(1 + 4*(trialID-1) : end) = h_POS(1 + 4*(trialID-1) : end) - trial_DUR(trialID);
            
        s_bool(start_POS : end_POS) = 1;
    end

    h_bool(1 + 4*(trialID-1) : 4*(trialID)) = trial_flag(trialID);
end

% Remove samples and trials contaminated by artifacts
s(logical(s_bool), :) = [];
h.EVENT.TYP(logical(h_bool), :) = [];
h.EVENT.DUR(logical(h_bool), :) = [];
h.EVENT.POS = h_POS;
h.EVENT.POS(logical(h_bool), :) = [];

end