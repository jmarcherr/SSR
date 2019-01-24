
fc = 1000; % carrier
fmod = 120; % fast modulation rate in hz
moddepth = 1; % of fast mod (for SAM)
stimrate = 2; % slow stimulus repetition rate in hz
fs = 48000;

Lstim = 2; % length of stimulus block
L = 3; % length of stimulus block incl pause
reps = 40; % how many block repetitions
Ltot = L*reps; % total length of stimulation experiment
fprintf('Total stimulation time: %d sec\n',Ltot)
stimtype = [1 0 0]; % AM tone, tone burst, chirp burst

t = 0:1/fs:Lstim-1/fs;

% create on/off slow modulation pattern:
%mod = max(square(2*pi*stimrate*t),0);
mod = zeros(1,Lstim*fs);
ons = 1:fs/stimrate:Lstim*fs+fs/stimrate;
Lr = .01*fs; % ramp length
r = [sin(linspace(0, pi/2, Lr)) ones(1,fs/stimrate/2-Lr*2) sin(linspace(pi/2, 0, Lr))]'; % ramp
for ii = 1:length(ons)-1
    mod(ons(ii):ons(ii)+fs/stimrate/2-1) = r;
end


% SAM tone
if stimtype(1)
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

% chirps
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

stim = stim.*mod; % apply slow on/off modulation 

stim(end:L*fs) = 0; % silence at the end of the block
t = 0:1/fs:L-1/fs;
plot(t,stim)
%soundsc(stim,fs)


%stim = repmat(stim,1,2);




