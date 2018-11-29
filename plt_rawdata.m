clear all
assr_startup
load(['JM2_data_demo2.mat'])


    cfgres = [];
    cfgres.resamplefs = 512;
    cfgres.detrend    = 'no';
    data_DG = ft_resampledata(cfgres, data_DG);
    
    data= data_DG;
 %%   
% [data_tmp,tcount,art_count] = thr_rjct(data,100,1,1,1);


 
%data = data_tmp;
%%
trigs = [1];% 2 3];
foihz = 1:1:100;
pow_cond={};
for k=1:length(trigs)
    %% Powspctrm
    idx=find(data.trialinfo==trigs(k))
    idx_tmp=idx;
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.foilim    = [foihz(1) foihz(end)];
    cfg.trials       = idx_tmp;
    cfg.keeptrials = 'no';
    cfg.channel      = {'Cz','FC6','FC5'};
    cfg.pad = 'maxperlen';
    cfg.output = 'pow';
    cfg.taper = 'dpss';
    cfg.tapsmofrq = 0.5;
    pow = ft_freqanalysis(cfg,data);
    fpow = [];
    fpow.label     = pow.label;
    fpow.freq      = pow.freq;
    fpow.dimord    = 'chan_freq_time';
    
    P = pow.powspctrm;   % copy the power spectrum
    fpow.P = P;
    
    pow_cond{k} = fpow;
end

%% Plotting power
close all
figure(1)
f = pow_cond{1,1}.freq;

for cc=1:length(pow_cond{1,1}.P(:,1))
for i=1:length(trigs)
    %subplot(3,1,i)
    loglog(f,squeeze(pow_cond{1,i}.P(cc,:))./(1./f),'linewidth',2)
    hold on
    xlim([3 100])
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
%%
%%%%%%%%%%%
%% Phase coherence (etire spectrum)

%foihz = [1 2 4 6 8 40];
trigs = [1];% 2 3];
foihz = 1:0.5:100;
chan_oi = [4 5 6 9 10    11    39    40    41    44    45    46  38 47];
for ff  = 1 : length(foihz)
    for k=1:length(trigs)
        %% ITPC
        idx=find(data.trialinfo==trigs(k))
        idx_tmp=idx;
        cfg = [];
        cfg.method = 'wavelet';
        cfg.toi    = -3:0.01:5;
        cfg.foi    = foihz(ff);
        cfg.trials       = idx_tmp;
        cfg.output = 'fourier';
        cfg.channel      = {'Cz','FC6','FC5'};
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
        itc.F = abs(F);
        
        % compute inter-trial phase coherence (itpc)
        itc.itpc      = F./abs(F);         % divide by amplitude
        itc.itpc      = sum(itc.itpc,1);   % sum angles
        itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
        itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
        
        % compute inter-trial linear coherence (itlc)
        itc.itlc      = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
        itc.itlc      = abs(itc.itlc);     % take the absolute value, i.e. ignore phase
        itc.itlc      = squeeze(itc.itlc); % remove the first singleton dimension
        itc_cond{k,ff} = itc;
        
    end
end
cd(rootdir)


%% Plotting section
% ITC
%Gather all
itc_all=[];
itc = [];
for ff=1:length(itc_cond)
    for kk=1:1
        itc(kk,ff,:,:,:) = itc_cond{kk,ff}.itpc;
        freqz(ff)= itc_cond{kk,ff}.freq;
    end
end

itc_all = squeeze(nanmean(itc,3)); % average over chan
itc_all(isnan(itc_all)) = 0;
close all

figure(2)
%titels = {'ITPC delta','ITPC gamma'}


cmp = cbrewer('seq','Blues',1000,'cubic');
for i=1:1
    subplot(3,1,i)
    h=imagesc(itc_cond{i}.time,freqz,itc_all(:,:))%,[0.2 0.5])
    axis xy
    %set(gca,'YScale','log');
    colormap(cmp)
    xlabel('time(s)')
    ylabel('Freq (Hz)')
    set(gca,'fontsize',16,'ytick',[4,40,205])
    xlim([-1 3.5])
    ylim([3 80])
    c=colorbar('location','EastOutside');
    ylabel(c, 'ITPC')
    set(gca,'fontsize',14)
    hold on
    plot(zeros(1,1000),linspace(0,300,1000),'k--')
    plot(3*ones(1,1000),linspace(0,300,1000),'k--')
end



foi = [4 40];
subplot(3,1,2)

p1=semilogx(freqz,nanmean(itc_all(:,:),2),'linewidth',2);
hold on
%p2=semilogx(freqz,squeeze(nanmean(itc_all(2,:,:),3)),'linewidth',2);
set(gca,'fontsize',14,'xtick',[0,4,10,20,30,40,50,60,205])
%xlabel('Freq(Hz)')
ylabel('ITPC')
ylim([0 0.6])
xlim([3 80])
plot(4*ones(1,1000),linspace(0,10,1000),'k--')
plot(40*ones(1,1000),linspace(0,10,1000),'k--')
%hleg=legend([p1],'4Hz delta','location','best')
hleg.Box = 'off'
grid on


set(gcf,'position',[481 112 464 693])

% power

subplot(3,1,3)

f = pow_cond{1,1}.freq;


for i=1:1
    loglog(f,squeeze(nanmean(pow_cond{1,i}.P,1))./(1./f),'linewidth',2)
    hold on
    xlim([3 80])
    ylim([0 100])
    %plot(4*ones(1,1000),linspace(0,10,1000),'k--')
    %plot(40*ones(1,1000),linspace(0,10,1000),'k--')
    set(gca,'fontsize',14,'xtick',[0,4,10,20,30,40,50,60,205])
    xlabel('Freq(Hz)')
    ylabel('Power (muV)')
    grid on
end
%toc
