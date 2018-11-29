close all
clear
rootdir = cd;
datadir = [rootdir,'/Speedmode_test'];
cd(datadir)
pow_cond={};
ids = [4 7];
for i=1:2
load(['40Hz_SP' num2str(ids(i)) '.mat']); %Minimal preprocessing
data = data_DG
%chan_oi = [7];
trigs = [1 2];
foihz = 1:1:80;
for kk=1:length(trigs)
    %% Powspctrm
    idx=find(data.trialinfo==trigs(kk))
    idx_tmp=idx;
    cfg = [];
    cfg.method = 'mtmfft';
    %cfg.toi    = -3:0.01:5;
    cfg.foilim    = [foihz(1) foihz(end)];
    cfg.trials       = idx_tmp;
    cfg.keeptrials = 'no';
    cfg.channel      = {'Fz', 'Cz'};
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
    
    pow_cond{i,kk} = fpow;
end
end
%%
close all
figure(1)
f = pow_cond{1}.freq;

for i=1:2
    subplot(2,1,i-1+kk-1)
    loglog(f,squeeze(pow_cond{i,1}.P)./(1./f),'linewidth',2)
    hold on
    loglog(f,squeeze(pow_cond{i,2}.P)./(1./f),'linewidth',2)
    xlim([3 500])
    ylim([-1 100])
    %plot(4*ones(1,1000),linspace(0,10,1000),'k--')
    %plot(40*ones(1,1000),linspace(0,10,1000),'k--')
    set(gca,'fontsize',14,'xtick',[0,4,10,20,30,40,50,60,205])
    xlabel('Freq(Hz)')
    ylabel('Power (muV)')
    grid on
    
end

% legend(data.label)
% cd(rootdir)


%% Phase coherence (etire spectrum)

%foihz = [1 2 4 6 8 40];
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
        cfg.channel      = {'Fz', 'Cz'};
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
    for kk=1:2
        itc(kk,ff,:,:,:) = itc_cond{kk,ff}.itpc;
        freqz(ff)= itc_cond{kk,ff}.freq;
    end
end

itc_all = squeeze(nanmean(itc,3)); % average over chan
close all

figure(2)
titels = {'ITPC delta','ITPC gamma'}

for kk=1:2
cmp = cbrewer('seq','Blues',1000,'cubic');
for i=1
    subplot(3,2,kk)
    h=contourf(itc_cond{i}.time,freqz,squeeze(itc_all(kk,:,:)),[0.25 0.6])
    axis xy
    set(gca,'YScale','log');
    colormap(cmp)
    xlabel('time(s)')
    ylabel('Freq (Hz)')
    set(gca,'fontsize',16,'ytick',[4,10,40,205])
    title(titels{i})
    xlim([-1 3.5])
    ylim([3 300])
    c=colorbar('location','EastOutside');
    %ylabel(c, 'ITPC')
    set(gca,'fontsize',14)
    hold on
    plot(zeros(1,1000),linspace(0,300,1000),'k--')
    plot(3*ones(1,1000),linspace(0,300,1000),'k--')
end




%%
foi = [4 40];
subplot(3,2,[3,4])

p1=semilogx(freqz,squeeze(nanmean(itc_all(kk,:,:),3)),'linewidth',2);
hold on
p2=semilogx(freqz,squeeze(nanmean(itc_all(kk,:,:),3)),'linewidth',2);
set(gca,'fontsize',14,'xtick',[0,4,10,20,30,40,50,60,205])
%xlabel('Freq(Hz)')
ylabel('ITPC')
ylim([0 0.6])
xlim([3 300])
plot(4*ones(1,1000),linspace(0,10,1000),'k--')
plot(40*ones(1,1000),linspace(0,10,1000),'k--')
%hleg=legend([p1 p2],'4Hz delta','40Hz gamma','location','best')
hleg.Box = 'off'
grid on
end