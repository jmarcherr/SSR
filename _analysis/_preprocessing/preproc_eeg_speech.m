clear all
ft_defaults
assr_startup

cd(bdfdir)
cd('./tech_aud/speech')
d = dir('*.bdf');
%load('./info/conditions.mat')
%% ------------Import data ----------------------------------------


for dd=1:length(d)
    cd(bdfdir)
    cd('./tech_aud/speech')
    dataset = d(dd).name;
    %% ------------Event extraction --------------------------------------
    
    triggers = [100];
    hdr = ft_read_header(dataset);
    cfg=[];
    cfg.layout =  'biosemi64.lay';
    cfg.continuous = 'yes';
    cfg.dataset = dataset;
    cfg.trialdef.eventtype    = 'STATUS';
    cfg.trialdef.eventvalue   = triggers;
    cfg.trialdef.prestim      = 0;
    cfg.trialdef.poststim     = 60;
    cfg = ft_definetrial(cfg);
    
    for tt=1:length(triggers)
        cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt;
    end
    
    
    %inital preprocessing(all channels no reref)
    data_int = ft_preprocessing(cfg);
    
    cd(datadir)
    cd('./tech_aud/speech')
    %Rereferencing (scalp)
    %cfg = [];
    cfg.dataset = dataset;
    cfg.channel     = {'Fp1','Fp2','Fz','Cz','T8','T7','FC5','FC6','EXG1','EXG2', '-Status'};%chaoi;
    cfg.reref       = 'yes';
    cfg.refchannel  = {'EXG1','EXG2'};
    cfg.layout      =  'biosemi64.lay';
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
    savefile = [dataset(1:end-4) '.mat']
    save(savefile,'data_DG','-v7.3');
    
    clear data_DG
    cd(rootdir)
end

