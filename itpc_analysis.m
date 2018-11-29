clear all
assr_startup
<<<<<<< HEAD
load(['JF_data_demo.mat'])

=======
addpath('Speedmode_test')
load(['40Hz_SP7.mat'])
addpath('C:\toolbox\fieldtrip_new')
addpath('C:\toolbox\cbrewer\')
ft_defaults
>>>>>>> d7f72d4fbb44eb50d26d5da7aa33dd5a49e3c1dd

    cfgres = [];
    cfgres.resamplefs = 512;
    cfgres.detrend    = 'no';
    data_DG = ft_resampledata(cfgres, data_DG);
    
    data= data_DG;
%%
trigs = [1];% 2 3];
<<<<<<< HEAD
foihz = 20:1:50;
=======
foihz = 1:1:100;
>>>>>>> d7f72d4fbb44eb50d26d5da7aa33dd5a49e3c1dd
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
    cfg.channel      = 'eeg';
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
<<<<<<< HEAD
    loglog(f,squeeze(pow_cond{1,i}.P(5,:)/pow_cond{1,i}.P(5,:))./(1./f),'linewidth',2)
=======
    loglog(f,squeeze(pow_cond{1,i}.P(:,:))./(1./f),'linewidth',2)
>>>>>>> d7f72d4fbb44eb50d26d5da7aa33dd5a49e3c1dd
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
%%
%%%%%%%%%%%
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
<<<<<<< HEAD
        cfg.channel      = chan_oi;
=======
        cfg.channel      = 'eeg';
>>>>>>> d7f72d4fbb44eb50d26d5da7aa33dd5a49e3c1dd
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
<<<<<<< HEAD
    for kk=1:2
=======
    for kk=1:1
>>>>>>> d7f72d4fbb44eb50d26d5da7aa33dd5a49e3c1dd
        itc(kk,ff,:,:,:) = itc_cond{kk,ff}.itpc;
        freqz(ff)= itc_cond{kk,ff}.freq;
    end
end

itc_all = squeeze(nanmean(itc,3)); % average over chan
close all

figure(2)
titels = {'ITPC delta','ITPC gamma'}


cmp = cbrewer('seq','Blues',1000,'cubic');
<<<<<<< HEAD
for i=1:2
    subplot(3,2,i)
    h=contourf(itc_cond{i}.time,freqz,squeeze(itc_all(i,:,:)),[0.25 0.4])
=======
for i=1
    subplot(3,1,i)
    h=contourf(itc_cond{i}.time,freqz,squeeze(itc_all(:,:)),[0 0.15])
>>>>>>> d7f72d4fbb44eb50d26d5da7aa33dd5a49e3c1dd
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

<<<<<<< HEAD
p1=semilogx(freqz,squeeze(nanmean(itc_all(1,:,:),3)),'linewidth',2);
hold on
p2=semilogx(freqz,squeeze(nanmean(itc_all(2,:,:),3)),'linewidth',2);
=======
p1=semilogx(freqz,squeeze(nanmean(itc_all(:,:),2)),'linewidth',2);
hold on
p2=semilogx(freqz,squeeze(nanmean(itc_all(:,:),2)),'linewidth',2);
>>>>>>> d7f72d4fbb44eb50d26d5da7aa33dd5a49e3c1dd
set(gca,'fontsize',14,'xtick',[0,4,10,20,30,40,50,60,205])
%xlabel('Freq(Hz)')
ylabel('ITPC')
ylim([0 0.6])
xlim([3 300])
plot(4*ones(1,1000),linspace(0,10,1000),'k--')
plot(40*ones(1,1000),linspace(0,10,1000),'k--')
<<<<<<< HEAD
hleg=legend([p1 p2],'4Hz delta','40Hz gamma','location','best')
=======
%hleg=legend([p1 p2],'4Hz delta','40Hz gamma','location','best')
>>>>>>> d7f72d4fbb44eb50d26d5da7aa33dd5a49e3c1dd
hleg.Box = 'off'
grid on


set(gcf,'position',[481 112 759 693])
%%
% power

subplot(3,2,[5,6])

f = pow_cond{1,1}.freq;


<<<<<<< HEAD
for i=1:2
=======
for i=1:1
>>>>>>> d7f72d4fbb44eb50d26d5da7aa33dd5a49e3c1dd
    loglog(f,squeeze(nanmean(pow_cond{1,i}.P,1))./(1./f),'linewidth',2)
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
<<<<<<< HEAD
toc
=======
%toc
>>>>>>> d7f72d4fbb44eb50d26d5da7aa33dd5a49e3c1dd
