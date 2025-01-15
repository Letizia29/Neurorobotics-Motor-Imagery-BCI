function [PSD, h_PSD, f] = get_PSD(s, h)

% This function computes the PSD given a signal s and its h file

fs = h.SampleRate;

% PSD
wlength = 0.5;              % [s]
pshift = 0.25;              % [s]
wshift = 0.0625;            % [s]
samplerate = fs;
mlength = 1;                % [s]

% PSD calculation
[PSD, f] = proc_spectrogram(s, wlength, wshift, pshift, samplerate, mlength);

% Select the frequencies from 4Hz to 48Hz
f_idx = (f >= 4 & f <= 48);
f = f(f_idx);

% Select the PSD from 4Hz to 48Hz
PSD = PSD(:, f_idx, :);

% Recompute POS and DUR

winconv = 'backward';
POS = proc_pos2win(h.EVENT.POS, wshift*h.SampleRate, winconv, wlength*h.SampleRate);

DUR = diff(POS) - 1;
DUR(length(POS)) = size(PSD, 1) - POS(length(POS)) - 1;

h_PSD.POS = POS;
h_PSD.DUR = DUR;
h_PSD.TYP = h.EVENT.TYP;

end