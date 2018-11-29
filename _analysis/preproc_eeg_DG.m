clear all
ft_defaults
assr_startup

    %% ------------Import data -----------------------------------------
    addpath([rootdir,'/_data']);
    dataset = ['JF_DGERP.bdf'];
    ft_defaults
    
%     %% ------------Bad channels (per subject)  --------------------------
%     clear badchans
%     badchans = badchans(s);

%chaoi = [32,11,20,7,24];
chaoi = {'Cz','P9','P10','T7','T8','FC5','FC6'}
    %% ------------Event extraction --------------------------------------
    triggers = [192,224];
    hdr = ft_read_header(dataset);
    cfg=[];
    cfg.layout =  'biosemi64.lay'; % why not 64?
    cfg.continuous = 'yes';
    cfg.dataset = dataset;
    cfg.trialdef.eventtype    = 'STATUS';
    cfg.trialdef.eventvalue   = triggers;
    cfg.trialdef.prestim      = 3;
    cfg.trialdef.poststim     = 5;
    cfg = ft_definetrial(cfg);
    
    for tt=1:length(triggers)
    cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt; 
    end

    
    %inital preprocessing(all channels no reref)
    cfg.channel = chaoi;
    data_int = ft_preprocessing(cfg);
    
    %Rereferencing (scalp)
    %cfg = [];
    cfg.dataset = dataset;
    cfg.channel     = chaoi;
    cfg.reref       = 'yes';
    cfg.refchannel  = {'T7','T8'}%setxor(data_int.label(1:64),badchans); % evt re-ref after channel removal
    cfg.layout      =  'biosemi64.lay';
    cfg.continuous  = 'yes';
    cfg.dftfilter   = 'yes';
    cfg.dftfreq     = [50, 100, 150, 200];
    cfg.hpfilter    = 'yes';
    cfg.hpfilttype  = 'firws';
    cfg.hpfreq      = 0.5;
    %cfg.lpfilter    = 'yes';
    %cfg.lpfilttype  = 'firws';
    %cfg.lpfreq      = 80;
    cfg.detrend     = 'yes';
    
    % rereferenced data struct
    data_DG = ft_preprocessing(cfg);
    
   
    % channel repair (ERP)
%     cfg = [];
%     cfg.layout        = 'biosemi64.lay';
%     cfg.method        = 'template';
%     cfg.channel       = data_int.label;
%     neighbours        = ft_prepare_neighbours(cfg);
%     cfg.method        = 'average';
%     cfg.layout        = 'biosemi64.lay';
%     cfg.neighbours    = neighbours;
% %    cfg.missingchannel   = badchans;
%     cfg.feedback      = 'no'
    
%     data_DG          = ft_channelrepair(cfg,data_DG);
%     cfg=[];
%     cfg.trials = 1;
%     tmp= ft_preprocessing(cfg,data_int)
%     
%     cfg =[];
%     tmp = ft_appenddata(cfg, tmp,data_DG)
%     
%     cfg=[];
%     cfg.trials = 2:length(tmp.trial);
%     tmp= ft_preprocessing(cfg,tmp)
%     % Check
%     
%     tmp.label
%     data_DG = tmp;
    
%     % Resample to 128Hz
%     cfgres = [];
%     cfgres.resamplefs = 128;
%     cfgres.detrend    = 'no';
%     data_DG = ft_resampledata(cfgres, data_DG);
    
    clear data_int tmp data_DG_ear
%    data_DG.missingchannels = badchans;
    % Data structure from this part:
    % data_erp (scalp data ERP + Ear-eeg)
    
    
    %%  Save mat
    
    save(['JF_data_demo.mat'],'data_DG','-v7.3')
    
    clear data_DG
    
    
    cd(rootdir)
%end

    % Resample to 1024Hz
%     cfgres = [];
%     cfgres.resamplefs = 1024;
%     cfgres.detrend    = 'no';
%     data_DG = ft_resampledata(cfgres, data_DG);
