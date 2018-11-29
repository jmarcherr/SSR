clear all
ft_defaults
assr_startup

    %% ------------Import data ----------------------------------------

    dataset = ['JM2_DG.bdf'];
    ft_defaults
    
    %% ------------Event extraction --------------------------------------
    cd(datadir)
    triggers = [150,120,160];
    hdr = ft_read_header(dataset);
    cfg=[];
    cfg.layout =  'biosemi32.lay'; % why not 64?
    cfg.continuous = 'yes';
    cfg.dataset = dataset;
    cfg.trialdef.eventtype    = 'STATUS';
    cfg.trialdef.eventvalue   = triggers;
    cfg.trialdef.prestim      = 0;
    cfg.trialdef.poststim     = 3;
    cfg = ft_definetrial(cfg);
    
    for tt=1:length(triggers)
    cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt; 
    end

    
    %inital preprocessing(all channels no reref)
    data_int = ft_preprocessing(cfg);
    
    
    %Rereferencing (scalp)
    %cfg = [];
    cfg.dataset = dataset;
    cfg.channel     = 'eeg';%chaoi;
    cfg.reref       = 'no';
    %cfg.refchannel  = {'T7','T8'}%setxor(data_int.label(1:64),badchans); % evt re-ref after channel removal
    cfg.layout      =  'biosemi32.lay';
    cfg.continuous  = 'yes';
    cfg.dftfilter   = 'yes';
    cfg.dftfreq     = [50, 100, 150, 200];
    cfg.hpfilter    = 'yes';
    cfg.hpfilttype  = 'firws';
    cfg.hpfreq      = 0.5;
    cfg.lpfilttype  = 'firws';
    
    % rereferenced data struct
    data_DG = ft_preprocessing(cfg,data_int);
    
    clear data_int 
    
    % Resample to 1024Hz
    cfgres = [];
    cfgres.resamplefs = 1024;
    cfgres.detrend    = 'no';
    data_int = ft_resampledata(cfgres, data_int);
    
    %%  Save mat
    save([dataset(1:end-4) '.mat'],'data_DG','-v7.3')    
    clear data_DG  
    cd(rootdir)

