clear all
ft_defaults
assr_startup

    %% ------------Import data ----------------------------------------
    cd(bdfdir)
    %dataset = ['APG_ASSR_IO_ModDepth_Ear-Right_Lvl-81dB_Mod-0.85.bdf'];
    dataset = 'pilot3_2xsamtone+oldx10min.bdf'
    ft_defaults
    
    %% ------------Event extraction --------------------------------------
    
    triggers = [100,200,255];
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
    data_int = ft_preprocessing(cfg);
    
    cd(datadir)
    %Rereferencing (scalp)
    %cfg = [];
    cfg.dataset = dataset;
    cfg.channel     = 'all';%chaoi;
    cfg.reref       = 'yes';
    cfg.refchannel  = {'EXG1','EXG2'}%setxor(data_int.label(1:64),badchans); % evt re-ref after channel removal
    cfg.layout      =  'biosemi32.lay';
    cfg.continuous  = 'yes';
    cfg.dftfilter   = 'yes';
    cfg.dftfreq     = [50, 100, 150, 200];
    cfg.lpfilttype  = 'but';
    cfg.lpfilter    = 'yes';
    cfg.lpfiltord   = 2;
    cfg.lpfreq      = 300;
    cfg.hpfilter    = 'yes';
    cfg.hpfilttype  = 'but'; 
    cfg.hpfreq      = 0.5;
    cfg.hpfiltord   = 2;

    
    % rereferenced data struct
    data_DG = ft_preprocessing(cfg,data_int);
    
    clear data_int 
    
    % Resample to 1024Hz
    cfgres = [];
    cfgres.resamplefs = 1024;
    cfgres.detrend    = 'no';
    data_DG = ft_resampledata(cfgres, data_DG);
    
    %%  Save mat
    savefile = [dataset(1:end-4) '_itpc.mat']
    save(savefile,'data_DG','-v7.3');    
    
    clear data_DG  
    cd(rootdir)

