clear all
assr_startup
ft_defaults

load(['pilot3_2xsamtone+oldx10min.mat']);

%%
%     cfgres = [];
%     cfgres.resamplefs = 512;
%     cfgres.detrend    = 'no';
%     data_DG = ft_resampledata(cfgres, data_DG);
%     
    data_all= data_DG;
    

%%
    idx=find(data_all.trialinfo==1)
    idx_tmp=idx;
    cfg = [];
    cfg.trials = idx_tmp;
    data_tmp = ft_selectdata(cfg,data_all);
    
    
%% own fft (concat)
f_fft = [];fs=data_tmp.fsample;
M =[];
%for tt = 1:1:length(data_tmp.trial)-10
M = cell2mat(data_tmp.trial(10:40));
plot(M(1,:)')
%pause(1.5)
%end
%n=nextpow2(length(M));
%f_fft = fft(M,2^n,2);
f_fft=fft(M)./(size(M,2)/2);

%%
%close all
figure(2)
for i=find(strcmp(data_tmp.label,'Cz'))
f = linspace(0,fs/2,length(f_fft)./2);
plot(f,abs(f_fft(i,1:end/2).^2))
hold on
end
%%
figure(3)
f = linspace(0,fs/2,length(f_fft)./2);
%plot(f,abs(nanmean(f_fft(find(strcmp(data_tmp.label,'Cz')),1:end/2),1)).^2)
plot(f,abs(f_fft(find(strcmp(data_tmp.label,'Cz')),1:end/2)).^2);
xlim([3 100])

%%

%%
trigs = [3];% 2 3];
%data_all.trialinfo = ones(size(data_all.trialinfo));
data = data_tmp;
foihz = 1:1:100;
pow_cond={};
for k=1:length(trigs)
    %% Powspctrm
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.foilim    = [foihz(1) foihz(end)];
    cfg.keeptrials = 'no';
    cfg.channel      = 'Cz';
    cfg.pad = 'nextpow2';
    cfg.output = 'pow';
    cfg.taper = 'hanning';
    cfg.tapsmofrq = .1;
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
    plot(f,squeeze(pow_cond{1,i}.P(:,:)),'linewidth',2)
    hold on
    xlim([3 60])
    %ylim([0 1])
    %plot(4*ones(1,1000),linspace(0,10,1000),'k--')
    %plot(40*ones(1,1000),linspace(0,10,1000),'k--')
    set(gca,'fontsize',14,'xtick',[0,4,10,20,30,40,50,60,205])
    xlabel('Freq(Hz)')
    ylabel('Power (muV)')
    grid on
end
end
legend(pow_cond{1}.label)
%%
%%%%%%%%%%%
%% Phase coherence (etire spectrum)

foihz = [1 2 4 6 8 40];
trigs=1;
%for ff  = 1 : length(foihz)
    
    for k=1%:length(trigs)
        %% ITPC
        idx=find(data_all.trialinfo==trigs(k))
        idx_tmp=idx;
        cfg = [];
        cfg.method = 'wavelet';
        cfg.toi    = 0:0.01:3;
        cfg.foi    = 4%foihz(ff);
        cfg.trials       = idx_tmp;
        cfg.output = 'fourier';
        cfg.channel      = 'all';
        cfg.pad = 'maxperlen';
        cfg.width = 12; %<- like dtu /hvidovre paper
        freq = ft_freqanalysis(cfg, data_all);
        
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
        %itc_cond{k,ff} = itc;
        
    end
%end
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
close all

figure(2)
titels = {'ITPC delta','ITPC gamma'}


cmp = cbrewer('seq','Blues',1000,'cubic');
for i=1
    subplot(3,1,i)
    h=contourf(itc_cond{i}.time,freqz,squeeze(itc_all(:,:)),[0 0.15])
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



foi = [4 40];
subplot(3,2,[3,4])

p1=semilogx(freqz,squeeze(nanmean(itc_all(:,:),2)),'linewidth',2);
hold on
p2=semilogx(freqz,squeeze(nanmean(itc_all(:,:),2)),'linewidth',2);
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


set(gcf,'position',[481 112 759 693])
%%
% power

subplot(3,2,[5,6])

f = pow_cond{1,1}.freq;


for i=1:1
    plot(f,squeeze(nanmean(pow_cond{1,i}.P,1)),'linewidth',2)
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
%toc
