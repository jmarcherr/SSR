function gamma_delta_ERP(audiogram)
%clear,clc

% Script to run gamma-rate (40 Hz) stimulation and ERP tone 

%
% Run ~5 min stimulation with:
% - stimrate = 40 ; L = 3 (gamma - 2 sec stimulation, 1 sec pause)

% Run with and without 4Hz delta bursts
% Run ~2 min stimulation with:
% - tone ERP stimrate 1Hz

% make sure that stimulus is presented to both ears!! (and not triggers)
deltas = [1 0 0]; %Delta yes or no
conds = [1 2 3];  %3 conditions (tone(delta/gamma) + ERP tone)
stimTrigger = [2 3 4]; %Stim trigger value for each condition
Ntrials = length(deltas); %Number of trials

gain_factor = 0.01;
%%%% compensation filter
origRMSlin = 10^(65/20);
adaptRMSlin = 10^((65+ (audiogram(3)/2))/20);
gain_factor = gain_factor.*(adaptRMSlin/origRMSlin);


for kk = 1:Ntrials
    % make stimulus:
    if kk<3 % if tone gamma/delta
        delta = deltas(kk);
        stimrate = 40; % in hz
        fs = 48000;
        L = 3; % length incl pause
        Lstim = 2; % length of stimulus block
        reps =60; % how many block repetitions
        Ltot = L*reps; % total length of stimulation experiment
        fprintf('Total stimulation time: %d sec\n',Ltot)
        Lr = round(.01*fs); % length of tone ramp
        
        tstim=0:1/fs:Lstim-1/fs;
        t=0:1/fs:L-1/fs;
        tone_f=1000;
        tone=(1+sin(stimrate*2*pi*tstim)).*sin(tone_f*2*pi*tstim);
                
        y_tone=zeros(L*fs,1);
        y_tone(1:Lstim*fs)=tone;
        ons = round(1:fs/stimrate:Lstim*fs); % click onsets in samples

        x(ons)=1;

        if length(unique(diff(ons)))>1
            warning('click rate does not match fs')
        end
        
        if delta
            ons4 = round(1:fs/4:Lstim*fs+fs/4);
            r=zeros(size(y_tone));
            for ii = 2:2:length(ons4)-1
                y_tone(ons4(ii):ons4(ii+1)) = 0;
                r (ons4(ii-1):ons4(ii))= [sin(linspace(0, pi/2, Lr)) ones(1,length(y_tone(ons4(ii):ons4(ii+1)))-Lr*2) sin(linspace(pi/2, 0, Lr))]';
                
            end
            y_tone = y_tone.*r;
        else

                    r = [sin(linspace(0, pi/2, Lr)) ones(1,Lstim*fs-Lr*2) sin(linspace(pi/2, 0, Lr))]';
                    y_tone(1:Lstim*fs) = y_tone(1:Lstim*fs).*r;
        end
        

            y=y_tone;

        %plot(y)

    %%
    % add trigger to y
    trig_dur = round(1e-3*fs); % trigger duration in samples
    admax = (2^8)-1;
    %stimTrigger = 2:(2); % triggers from 2-8 are useable
    trig_word = (2.^stimTrigger(kk))-1;
    trig_ampl = trig_word./admax;
    trigger = zeros(L*fs,1);
    trigger(1:trig_dur+1,1) = trig_ampl(1); % start trigger each block
    
    y=y.*gain_factor;
    y = [y y trigger]';
    y = repmat(y,1,reps);
    else
        % stimulus parameters for ERP tone
        stimrate = 1; % average tone rate, in hz (check this)
        fs = 48000;
        Lstim = 10; % length of stimulation
        toneL = .1; % length of tones
        toneF = 1000; % tone freq
        reps =18; % how many block repetitions
        Ltot = Lstim*reps; % total length of stimulation experiment
        jitmax = .25; % max temporal jitter, in sec
        Lr = round(.01*fs); % length of tone ramp
        fprintf('Total stimulation time: %d sec\n',Ltot)
        
        % generate one tone
        r = [sin(linspace(0, pi/2, Lr)) ones(1,fs*toneL-Lr*2) sin(linspace(pi/2, 0, Lr))];
        a = r.*sin(2*pi*toneF*linspace(0,toneL,toneL*fs));
        
        % find onsets
        ons = round(1:fs/stimrate:Ltot*fs); % in samples
        for oo = 2:length(ons)
            ons(oo) = ons(oo) +  round(fs*(2*jitmax*rand(1)-jitmax));
        end
        
        % make trigger
        % add trigger to y
        trig_dur = round(1e-3*fs); % trigger duration in samples
        admax = (2^8)-1;
        %stimTrigger = 2:(2); % triggers from 2-8 are useable
        trig_word = (2.^stimTrigger(kk))-1;
        trig_ampl = trig_word./admax;
        
        % make tone and trigger sequence
        y = zeros(Ltot*fs,1);
        trigger = zeros(Ltot*fs,1);
        for ii = 1:length(ons)
            y(ons(ii):ons(ii)+length(a)-1,:) = a;
            trigger(ons(ii):ons(ii)+trig_dur-1,1) = trig_ampl(1);
        end
        y=y.*gain_factor;
        y=[y y trigger]';
    end    
    
%% start psychportaudio    
     %init psychportaudio
    try PsychPortAudio('GetOpenDeviceCount')
        PsychPortAudio('close');
    end
    InitializePsychSound;
    dev = PsychPortAudio('GetDevices');
    if kk==1
        devid =6;
    else
        devid =0;
    end
    %devid  3; % 0 for PHYS2 HEAAUD
    selectchannel = [4 5 12 ;0 0 0]; %  [4 12;0 0]; ER2 + adat3
    nchans =  size(selectchannel,2);
    pah = PsychPortAudio('Open', devid, [], 1, fs, nchans, [], [], selectchannel,4);
    PsychPortAudio('FillBuffer', pah, y);%y(1:2,:));
    PsychPortAudio('Volume',pah,0);PsychPortAudio('Start', pah, 1, 0, 0, .2);PsychPortAudio('Stop', pah, 1);PsychPortAudio('Volume',pah,1);
    
    
    %% setup screen and keyboard
    if kk==1
        Screen('Preference', 'Verbosity', 0); % Suppress warnings from PTB3
        KbName('UnifyKeyNames');
        escKey = KbName('ESCAPE');
        quitKey = KbName('q');
        no_screens=Screen('Screens'); % external screens?
        thescreen =max(no_screens); % if so, choose external
        min_factor = 0.999999999999999999; % how much the screen should shrink
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
        
    end
    
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
        if kk<3
            DrawFormattedText(w, double(num2str('Klar til næste test? \n \n Tryk for at starte')), 'center','center');
            Screen('Flip',w,0,0,1);
            while KbCheck; end
            WaitSecs(.1);
            KbWait(-1);
        else
            DrawFormattedText(w, double(num2str('Forsøg slut')), 'center','center');
            Screen('Flip',w,0,0,1);
            while KbCheck; end
            WaitSecs(.1);
            KbWait(-1);
        end
    catch
        ShowCursor;
        ListenChar(0);
        Screen('CloseAll')
        psychrethrow(psychlasterror);
        PsychPortAudio('Stop', pah);
    end
end
ShowCursor;
ListenChar(0);
PsychPortAudio('Stop', pah);
Screen('CloseAll')


end