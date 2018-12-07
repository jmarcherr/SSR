clear;close all
addpath('_func')

%%
%Script to run ASSR experiment
%The following script can present different parameters of the ASSR stimulus
%designed to elicit a EFR-response at both low and high rates

gain_factor_calib = 0.01;
% corresponded to 77.8 LEeq (dB) 118.4 peSPL (13/11/18)-> for 70dB 
gain_factor_75 = gain_factor_calib*10^((75-77.8)/20);
% for 85dB
gain_factor_85 = gain_factor_calib*10^((85-77.8)/20);

fc = [2005]; % carrier
fmod =[93];% 205]; % fast modulation rate in hz
moddepth = [1]; % of fast mod (for SAM)
stimrate = [2]; % slow stimulus repetition rate in hz

fs = 48000;

Lstim = 2; % length of stimulus block
L = 3; % length of stimulus block incl pause
reps = 1;%0; % how many block repetitions
Ltot = L*reps; % total length of stimulation experiment
fprintf('Total stimulation time: %d sec\n',Ltot)
stimtype = [1 0 0]; % SAM tone, tone burst, chirp burst
gain_factor = gain_factor_75;
%%%% Make stimulus
y=[];
[y,cfg] = make_assr_stim(fc,fmod,moddepth,stimrate,fs,Lstim,L,stimtype,gain_factor);

%%%% add trigger to stim
trigger = [ones(1,round(1e-3 * fs)) zeros(1,length(y)-round(1e-3 * fs))]';

y=y.*gain_factor;
%%
y = [y y]%trigger]';
t=0:1/fs:length(y)/fs-1/fs;
plot(t,y)
%%
y = repmat(y,1,reps)';

%% start psychportaudio
%init psychportaudio
try PsychPortAudio('GetOpenDeviceCount')
    PsychPortAudio('close');
end
InitializePsychSound;
dev = PsychPortAudio('GetDevices');
devid = 1;
%devid  3; % 0 for PHYS2 HEAAUD
selectchannel = [1 2;0 0];% 12 ;0 0 0]; %  [4 12;0 0]; ER2 + adat3
nchans =  size(selectchannel,2);
pah = PsychPortAudio('Open', devid, [], 1, fs, nchans, [], [], selectchannel,4);
PsychPortAudio('FillBuffer', pah, y');%y(1:2,:));
PsychPortAudio('Volume',pah,0);PsychPortAudio('Start', pah, 1, 0, 0, .2);PsychPortAudio('Stop', pah, 1);PsychPortAudio('Volume',pah,1);

%% initialize trigger box communication

% trig = HEATriggerbox();
% trig.find_triggerbox_win();
% trig.connect();
% if trig.is_connected()
%     trig.set_trigger(10);
% end
%% setup screen and keyboard

Screen('Preference', 'Verbosity', 0); % Suppress warnings from PTB3
KbName('UnifyKeyNames');
escKey = KbName('ESCAPE');
quitKey = KbName('q');
no_screens=Screen('Screens'); % external screens?
thescreen =max(no_screens); % if so, choose external
min_factor =0.999999999999999999; % how much the screen should shrink
scrndims = Screen(thescreen,'rect')*min_factor;
full_scrn = Screen(thescreen,'rect');
scrndims = [4 5 scrndims(3) scrndims(4)]%[full_scrn(3)-scrndims(3) full_scrn(4)-scrndims(4) scrndims(3) scrndims(4)];


rect = CenterRect([0 0 scrndims(3)-(full_scrn(3)-scrndims(3)) scrndims(4)-(full_scrn(4)-scrndims(4))],scrndims);
backgroundcolor = [100 100 100];
black = [0 0 0];
white = [255 255 255];
HideCursor;
%Experiment window
% open window
[w,rect]=Screen('OpenWindow', thescreen, backgroundcolor, rect); % open main window
Screen('TextSize',w, 40);
%Screen('TextFont',w, 'New Courier Bold');
Screen('TextColor',w,[255 255 255]); % white
HideCursor;

% ready window
DrawFormattedText(w, double(num2str('Forsøg klar \n \n Tryk for at starte')), 'center','center');
Screen('Flip',w,0,0,1);
while KbCheck; end
WaitSecs(.1);
KbWait(-1);


DrawFormattedText(w, double(num2str('+')), 'center','center');
Screen('Flip',w,0,0,1);


 %% RUN STIMULATION

KbName('UnifyKeyNames');
escKey = KbName('ESCAPE');
commandwindow;
try
    
    % ListenChar(2);
    
    breakflag = 0;
    offset = 1;
    fprintf('Begin stimulation\n')
    fprintf('Press ESC to exit\n')
    tnow = GetSecs;
    PsychPortAudio('Start', pah, 1, tnow+offset);
    
    s = PsychPortAudio('GetStatus', pah);
    if ~s.Active
        while 1 % wait untill audio starts
            s = PsychPortAudio('GetStatus', pah);
            [~,~,keyCode] = KbCheck(-1);
            if s.Active
                break;
            end
            if keyCode(escKey)
                fprintf('escape\n')
                breakflag=1;
                PsychPortAudio('Stop', pah, 0);
                break;
            end
        end
    end
    
    s = PsychPortAudio('GetStatus', pah);
    if s.Active
        tStart = GetSecs;
        while 1 % wait untill audio stops
            s = PsychPortAudio('GetStatus', pah);
            [~,~,keyCode] = KbCheck(-1);
            if ~s.Active
                break
            end
            if keyCode(escKey)
                fprintf('Escape\n')
                breakflag=1;
                PsychPortAudio('Stop', pah, 0);
                break
            end
        end
    end
    
    PsychPortAudio('Stop', pah, 1);
    
    fprintf('End of stimulation\n')
    
    
    DrawFormattedText(w, double(num2str('Forsøg slut')), 'center','center');
    Screen('Flip',w,0,0,1);
    while KbCheck; end
    WaitSecs(.1);
    KbWait(-1);
    
catch
    ShowCursor;
    ListenChar(0);
    Screen('CloseAll')
    psychrethrow(psychlasterror);
    PsychPortAudio('Stop', pah);
end

ShowCursor;
ListenChar(0);
PsychPortAudio('Stop', pah);
Screen('CloseAll')


%end