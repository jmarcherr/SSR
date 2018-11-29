function [stim,cfg] = make_assr_stim(fc,fmod,moddepth,stimrate,fs,Lstim,L,stimtype,gain_factor)
%make_assr_stim: Function to generate steady-state stimulus with slow, and
%fast modulation rates.
%   fc          : Carrier frequency (Hz)
%   fmod        : Modulation frequency (Hz)
%   moddepth    : Modulation depth (%)
%   stimrate    : Slow modulation frequency (s)
%   fs          : Sampling frequency
%   Lstim       : Length of stimulation period
%   L           : Length of trial
%   stimtype    : type of stimulus (1)SAM-tone (2)Tones (3)Chirps

cfg = struct;
cfg.fc = fc;
cfg.fmod = fmod;
cfg.moddepth = moddepth;
cfg.stimrate = stimrate;
cfg.fs = fs;
cfg.Lstim = Lstim;
cfg.L = L;
cfg.stimtype = stimtype;
cfg.gain_factor = gain_factor;
t = 0:1/fs:Lstim-1/fs;


% create on/off slow modulation pattern:
mod = zeros(1,Lstim*fs);
ons = 1:fs/stimrate:Lstim*fs+fs/stimrate;
Lr = .01*fs; % ramp length
r = [sin(linspace(0, pi/2, Lr)) ones(1,fs/stimrate/2-Lr*2) sin(linspace(pi/2, 0, Lr))]'; % ramp
for ii = 1:length(ons)-1
    mod(ons(ii):ons(ii)+fs/stimrate/2-1) = r;
end

% SAM tone
if stimtype(1)
    % Test for ERB modulation
    ERB_cb = (24.7*(4.37*(fc/1000)+1))/2;
    if fmod>ERB_cb
        warning('modulation rate too high for the chosen carrier frequency')
    end
    stim = (( 1 + moddepth*sin(2*pi*fmod*t))/2) .* sin(2*pi*fc*t);
end

% tones
if stimtype(2)
    stim = sin(2*pi*fc*t);
    tonemod = zeros(1,Lstim*fs);
    ons = 1:fs/fmod:Lstim*fs+fs/fmod;
    Lr = .001*fs; % ramp length
    r = [sin(linspace(0, pi/2, Lr)) ones(1,fs/fmod/2-Lr*2) sin(linspace(pi/2, 0, Lr))]'; % ramp
    for ii = 1:length(ons)-1
        tonemod(ons(ii):ons(ii)+fs/fmod/2-1) = r;
    end
    stim = stim.*tonemod;
end


% Double sam tone
if stimtype(3)
    stim = zeros(1,Lstim*fs);
    a = bmchirp(100,10000,fs,0,0);
    if fmod>fs/length(a)
        error('chirp too short for modulation rate')
    end
    ons = 1:fs/fmod:Lstim*fs+fs/fmod;
    for ii = 1:length(ons)-1
        stim(ons(ii):ons(ii)+length(a)-1) = a;
    end
end


% make trial
stim = stim.*mod; % apply slow on/off modulation
stim(end:L*fs) = 0; % silence at the end of the block
t = 0:1/fs:L-1/fs;

%plot(t,stim)
stim = stim';


end

