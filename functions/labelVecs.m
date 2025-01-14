function vectors = labelVecs(s, h)

% This function calculates the label vectors for a given signal s and
% its information structure h.

% Creation of the label vectors for the provided GDF file
% Create Tk [trial vector] (1, 2, 3, ... N)
% Create Fk [fixation periods] (0 or event value)
% Create Ak [cue periods] (0 or event value)
% Create CFk [continuous feedback periods] (0 or event value)
% Create Xk [hit/miss periods] (0 or event value)
% Create Wk [windows] for PSD

% events codes
trial_start = 1;
fixation = 786;
cue_hand = 773;
cue_feet = 771;
cont_feedback = 781;
hit = 897;
miss = 898;



% Initialize the data structures to store the information.
Tk = zeros(size(s, 1), 1); % Trial vector
Fk = zeros(size(s, 1), 1);
Ak = zeros(size(s, 1), 1);
CFk = zeros(size(s, 1), 1);
Xk = zeros(size(s, 1), 1);
Wk = zeros(size(s, 1), 1);

% computation
trial_number = 0;

for i = 1:length(h.TYP)
    % extracting trials from fixation to end of continuous feedback period
    if h.TYP(i,1) == fixation
        trial_number = trial_number + 1;
        start = h.POS(i,1);
        finish = start + h.DUR(i,1) + h.DUR(i+1,1) + h.DUR(i+2,1) - 1;
        Tk(start:finish,1) = trial_number;
    end 
   
    % fixation periods
    if h.TYP(i,1) == fixation
        start = h.POS(i,1);
        finish = start + h.DUR(i,1);
        Fk(start:finish,1) = fixation;
    end

    % cue periods
    if h.TYP(i,1) == cue_hand
        start = h.POS(i,1);
        finish = start + h.DUR(i,1);
        Ak(start:finish,1) = cue_hand;
    end
    if  h.TYP(i,1) == cue_feet  
        start = h.POS(i,1);
        finish = start + h.DUR(i,1);
        Ak(start:finish,1) = cue_feet;
    end
    
    % continuous feedback periods
    if h.TYP(i,1) == cont_feedback
        start = h.POS(i,1);
        finish = start + h.DUR(i,1);
        CFk(start:finish,1) = cont_feedback;
    end

    % hit/miss periods
    if h.TYP(i,1) == hit
        start = h.POS(i,1);
        finish = start + h.DUR(i,1);
        Xk(start:finish,1) = hit;
    end
    
    if h.TYP(i,1) == miss
        start = h.POS(i,1);
        finish = start + h.DUR(i,1);
        Xk(start:finish,1) = miss;
    end

    % extracting windows from cue to end of continuous feedback period
    if h.TYP(i,1) == cue_hand || h.TYP(i,1) == cue_feet
        start = h.POS(i,1);
        finish = start + h.DUR(i,1) + h.DUR(i+1,1) - 1;
        Wk(start:finish,1) = 1;
    end


    vectors = struct();
    vectors.Tk = Tk;
    vectors.Fk = Fk;
    vectors.Ak = Ak;
    vectors.CFk = CFk;
    vectors.Xk = Xk;
    vectors.Wk = Wk;


end