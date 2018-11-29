%% plot for subject JF in WP4 data
clear
rootdir = cd;
datadir = [rootdir,'/_data'];
cd(datadir)

load(['JF_data_demo.mat']); %Minimal preprocessing
data = data_DG
chan_oi = [7];
trigs = [1 2];
foihz = 1:1:300;
pow_cond={};
for k=1:2
    %% Powspctrm
    idx=find(data.trialinfo==trigs(k))
    idx_tmp=idx;
    cfg = [];
    cfg.method = 'mtmfft';
    %cfg.toi    = -3:0.01:5;
    cfg.foilim    = [foihz(1) foihz(end)];
    cfg.trials       = idx_tmp;
    cfg.keeptrials = 'no';
    cfg.output = 'fourier';
    cfg.channel      = 'eeg';
    cfg.pad = 'maxperlen';
    cfg.output = 'pow';
    cfg.taper = 'hanning';
    pow = ft_freqanalysis(cfg,data);
    fpow = [];
    fpow.label     = pow.label;
    fpow.freq      = pow.freq;
    fpow.dimord    = 'chan_freq_time';
    
    P = pow.powspctrm;   % copy the Fourier spectrum
    %N = size(F,1);           % number of trials
    
    fpow.P = P;
    
    pow_cond{k} = fpow;
end
close all
figure(1)
f = pow_cond{1,1}.freq;

for cc=1:length(pow_cond{1,1}.P(:,1))
for i=1:2
    subplot(3,1,i)
    loglog(f,squeeze(pow_cond{1,i}.P(cc,:))./(1./f),'linewidth',2)
    hold on
    xlim([3 300])
    ylim([0 100])
    %plot(4*ones(1,1000),linspace(0,10,1000),'k--')
    %plot(40*ones(1,1000),linspace(0,10,1000),'k--')
    set(gca,'fontsize',14,'xtick',[0,4,10,20,30,40,50,60,205])
    xlabel('Freq(Hz)')
    ylabel('Power (muV)')
    grid on
end
end
legend(data.label)
cd(rootdir)
