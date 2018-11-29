%% Test make_assr_stim
clear all
%parameters
fc = 1500; % carrier
fmod =80; % fast modulation rate in hz
moddepth = 1; % of fast mod (for SAM)
stimrate = 2; % slow stimulus repetition rate in hz
fs = 48000;

Lstim = 2; % length of stimulus block
L = 3; % length of stimulus block incl pause
reps = 40; % how many block repetitions
Ltot = L*reps; % total length of stimulation experiment
fprintf('Total stimulation time: %d sec\n',Ltot)
stimtype = [1 0 0]; % SAM tone, tone burst, chirp burst

[stim,cfg] = make_assr_stim(fc,fmod,moddepth,stimrate,fs,Lstim,L,stimtype);

%soundsc(stim,fs)