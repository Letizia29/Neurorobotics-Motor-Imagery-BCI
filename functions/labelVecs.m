function vectors = labelVecs(s, h)

% This function calculates the label vectors for a given signal s and
% its information structure h.

% Initialize the data structures to store the information.
Tk = zeros(1, size(s, 1)); % Trial vector
Fk = zeros(1, size(s, 1));
Ck = zeros(1, size(s, 1));
CFk = zeros(1, size(s, 1));
Xk = zeros(1, size(s, 1));

% Find when a trial begins or ends (1 for offline trial beginning, 897 or 898 for online data ending)
idx_type_init = ((h.EVENT.TYP == 1) | (h.EVENT.TYP == 897 | h.EVENT.TYP == 898));
num_trials = sum(idx_type_init); % Find the number of trials

% Identify the index of different events
idx_type_fix = find(h.EVENT.TYP == 786);
idx_type_cue = find(h.EVENT.TYP == 771 | h.EVENT.TYP == 773 | h.EVENT.TYP == 783);
idx_type_contfb = find(h.EVENT.TYP == 781);
idx_type_hitmiss = find(h.EVENT.TYP == 897 | h.EVENT.TYP == 898);

% Ifentify when a precise event happens and its duration
Tk_pos = h.EVENT.POS(idx_type_init);
Tk_dur = h.EVENT.DUR(idx_type_init);

Fk_pos = h.EVENT.POS(idx_type_fix);
Fk_dur = h.EVENT.DUR(idx_type_fix);

Ck_pos = h.EVENT.POS(idx_type_cue);
Ck_dur = h.EVENT.DUR(idx_type_cue);

CFk_pos = h.EVENT.POS(idx_type_contfb);
CFk_dur = h.EVENT.DUR(idx_type_contfb);

Xk_pos = h.EVENT.POS(idx_type_hitmiss);
Xk_dur = h.EVENT.DUR(idx_type_hitmiss);

% Build Tk
for i = 1 : num_trials-1
    Tk(Tk_pos(i) : Tk_pos(i+1)) = i;
end

Tk(Tk_pos(num_trials) : end) = num_trials;
% I did it here since in the for loop after this one the indeces exceeded
% the array size.

% Build Fk, Ck, CFk, Xk:
% The value of an element in each vector will be equal to the code of the
% specific event

pos_in_trial_Ck = find(h.EVENT.TYP == 773 | h.EVENT.TYP == 771, 1, 'first');
pos_in_trial_CFk = find(h.EVENT.TYP == 781, 1, 'first');
pos_in_trial_Xk = find(h.EVENT.TYP == 897 | h.EVENT.TYP == 898, 1, 'first');

% NB THERE ARE SOME PROBLEMS IN THE ONLINE DATA
for i = 1 : num_trials
    Fk(Fk_pos(i) : Fk_pos(i)+Fk_dur(i)-1) = 1;
    Ck(Ck_pos(i) : Ck_pos(i)+Ck_dur(i)-1) = h.EVENT.TYP(pos_in_trial_Ck+4*(i-1));       % Da migliorare ma funziona
    CFk(CFk_pos(i) : CFk_pos(i)+CFk_dur(i)-1) = h.EVENT.TYP(pos_in_trial_CFk+4*(i-1));

    if(isempty(idx_type_hitmiss) == false)          % Online file: the last event in a trial is a hit or miss
        Xk(Xk_pos(i) : Xk_pos(i)+Xk_dur(i)-1) = h.EVENT.TYP(pos_in_trial_Xk+4*(i-1));
    end
end


vectors.Tk = Tk;
vectors.Fk = Fk;
vectors.Ck = Ck;
vectors.CFk = CFk;
vectors.Xk = Xk;

end