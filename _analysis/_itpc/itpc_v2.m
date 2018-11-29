clear all
assr_startup

%for i=1:length(SUBJECTS)
i=1;
    s=SUBJECTS{i}
    datadir = ['_data'];
    cd(datadir)
    
    load([s '_data_DG.mat'])
    data= data_DG;
    
%% Denoising     
    % old: lp 20 Hz, 0-2 s latency; thr_rjct(data2,150,0.01,1,0);
    cfg = [];
    cfg.lpfreq = 30;
    cfg.lpfilter = 'yes';
    
    [data2] = ft_preprocessing(cfg, data)
    
    cfg = [];
    cfg.latency = [0 2];
    [data2] = ft_selectdata(cfg, data2);
    
   % if ear_eeg
    %    keep_trials = 1:length(data.trial);
    %else
        [data_clean,t_count] = thr_rjct(data2,150,0.01,1,0);
        
        
        keep_trials = [];
        for ti = 1 : length(data_clean.trial)
            keep_trials(ti)= length(find(isnan(data_clean.trial{ti}(:))))==0;
        end
        disp(mean(keep_trials))
        
        %
        cfg = [];
        cfg.trials = find(keep_trials)
        [data] = ft_selectdata(cfg, data);
   % end
    
%%    
    
    itc_cond = {};
    foihz = 1:0.5:60;
    %foihz = [1 2 4 6 8 40];
    for ff  = 1 : length(foihz)
        % Phase coherence
        trigs = [1 2];
        for k=1:length(trigs)
            clear freq itc
            idx=find(data.trialinfo==trigs(k));
            idx_tmp=idx;
            cfg = [];
            cfg.method = 'wavelet';
            cfg.toi    = -3:0.01:5;
            cfg.foi    = foihz(ff);
            cfg.trials       = idx_tmp;
            cfg.output = 'fourier';
            cfg.channel      = 'EEG';
            cfg.pad = 'maxperlen';
            cfg.width = 12; %<- like dtu /hvidovre paper
            freq = ft_freqanalysis(cfg, data);
            
            
            itc = [];
            itc.label     = freq.label;
            itc.freq      = freq.freq;
            itc.time      = freq.time;
            itc.dimord    = 'chan_freq_time';
            
            F = freq.fourierspctrm;   % copy the Fourier spectrum
            N = size(F,1);           % number of trials
            
            % compute inter-trial phase coherence (itpc)
            itc.itpc      = F./abs(F);         % divide by amplitude
            itc.itpc      = nansum(itc.itpc,1);   % sum angles
            itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
            itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
            
            %itc.itpowc      = abs(F);         % divide by amplitude
            %itc.itpowc      = sum(itc.itpowc,1);   % sum angles
            %itc.itpowc      = abs(itc.itpowc)/N;   % take the absolute value and normalize
            %itc.itpowc      = squeeze(itc.itpowc); % remove the first singleton dimension
            
            % compute inter-trial linear coherence (itlc)
            %itc.itlc      = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
            %itc.itlc      = abs(itc.itlc);     % take the absolute value, i.e. ignore phase
            %itc.itlc      = squeeze(itc.itlc); % remove the first singleton dimension
            itc.keep_trials = keep_trials;
            itc_cond{k,ff} = itc;
            
        end
        cd(rootdir)
        
    end
    % cd(fullfile(rootdir,'_itpc','results'))

        save([fullfile(datadir,SUBJECTS{i}),'_itpc.mat'],'itc_cond')

%end
%%

%%

%for i=1:length(SUBJECTS)
i=1
    s=SUBJECTS{i}
    datadir = [cd filesep s filesep 'Experiment'];
    cd(datadir)
    
    load([s '_data_DG2.mat'])
    data= data_DG;
    
    %[data,todss] = eog_denoising(data,1:64,0.8);
    
    
    itc_cond = {};
    foihz = [4 40]
    for ff  = 1 : length(foihz)
        % Phase coherence
        trigs = [1 2];
        for k=1:length(trigs)
            clear freq itc
            idx=find(data.trialinfo==trigs(k))
            idx_tmp=idx;
            cfg = [];
            cfg.method = 'wavelet';
            cfg.toi    = -1:0.01:3;
            cfg.foi    = foihz(ff);
            cfg.trials       = idx_tmp;
            cfg.output = 'fourier';
            cfg.channel      = 'EEG';
            cfg.pad = 'maxperlen';
            freq = ft_freqanalysis(cfg, data);
            
            
            itc = [];
            itc.label     = freq.label;
            itc.freq      = freq.freq;
            itc.time      = freq.time;
            itc.dimord    = 'chan_freq_time';
            
            F = freq.fourierspctrm;   % copy the Fourier spectrum
            N = size(F,1);           % number of trials
            
            % compute inter-trial phase coherence (itpc)
            itc.itpc      = F./abs(F);         % divide by amplitude
            itc.itpc      = sum(itc.itpc,1);   % sum angles
            itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
            itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
            itc_cond{k,ff} = itc;
            
        end
        cd(rootdir)
        
    end
    % cd(fullfile(rootdir,'_itpc','results'))
    
    save([fullfile(datadir,SUBJECTS{i}),'_itpc.mat'],'itc_cond')
%end