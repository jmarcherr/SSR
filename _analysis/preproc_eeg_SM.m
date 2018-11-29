clear all
ft_defaults
assr_startup
addpath('Speedmode_test')
    %% ------------Import data ----------------------------------------
%for i=[4 7]
    dataset = ['JM2_DG.bdf'];
    ft_defaults
    
%     %% ------------Bad channels (per subject)  --------------------------
%     clear badchans
%     badchans = badchans(s);

%chaoi = [32,11,20,7,24];
%chaoi = {'Cz','P9','P10','T7','T8','FC5','FC6'}
    %% ------------Event extraction --------------------------------------
    cd('_data')
    triggers = [150,120,160];
    hdr = ft_read_header(dataset);
    cfg=[];
    cfg.layout =  'biosemi32.lay'; % why not 64?
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
    data_int = ft_preprocessing(cfg);
        %Resample to 128Hz
    cfgres = [];
    cfgres.resamplefs = 512;
    data_int = ft_resampledata(cfgres, data_int);
    
    %Rereferencing (scalp)
    %cfg = [];
    cfg.dataset = dataset;
    cfg.channel     = 'eeg';%chaoi;
    cfg.reref       = 'no';
    cfg.refchannel  = {'T7','T8'};%setxor(data_int.label(1:64),badchans); % evt re-ref after channel removal
    cfg.layout      =  'biosemi32.lay';
    cfg.continuous  = 'yes';
    cfg.hpfilter    = 'yes';
    cfg.hpfilttype  = 'firws';
    cfg.hpfreq      = 0.5;
    %cfg.lpfreq      = 100;
    %cfg.lpfilter    = 'yes';
    %cfg.lpfilttype  = 'firws';
    %cfg.detrend     = 'yes';
    cfg.detrend    = 'no';
    cfg.dftfilter   = 'yes';
    cfg.dftfreq     = [50, 100, 150, 200];
    
    % rereferenced data struct
    data_DG = ft_preprocessing(cfg,data_int);
    

    
    clear data_int tmp data_DG_ear

    
    
    %%  Save mat
    
    save(['JM2_data_DG.mat'],'data_DG','-v7.3')
    
    clear data_DG
    
    
    cd(rootdir)
%end
