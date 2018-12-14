clear all
assr_startup
ft_defaults

load(['pilot3_2xsamtone+oldx10min_itpc.mat']);
data_all= data_DG;

%%%%%%%%%%%
%% Phase coherence (etire spectrum)

foihz = 100:300;%[1:7; 204:210];
trigs=[1,2,3];
for ff  = 1 : size(foihz,1)
    
    for k=1:length(trigs)
        %% ITPC
        idx=find(data_all.trialinfo==trigs(k))
        idx_tmp=idx;
        cfg = [];
        cfg.method = 'wavelet';
        cfg.toi    = -3:0.01:5;
        cfg.foi    = foihz(ff,:);
        cfg.trials       = idx_tmp;
        cfg.output = 'fourier';
        cfg.channel      = {'all', '-Status'};
        cfg.pad = 'maxperlen';
        %cfg.width = 12; %<- like dtu /hvidovre paper
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
        itc_cond{k,ff} = itc;
        
    end
end
cd(rootdir)

%% Plotting section
% ITC
%Gather all
itc_all=[];
itc = [];
for ff=1:1
    for kk=1:3
        itc(kk,ff,:,:,:) = itc_cond{kk,ff}.itpc;
        freqs(ff,:)= itc_cond{kk,ff}.freq;
    end
end

%itc_all = squeeze(itc,3); % average over chan
close all

figure(2)
titels = {'ITPC delta','ITPC gamma'}

ic=0;
cmp = cbrewer('seq','Blues',11,'cubic');
for ff=1%:2
for i=1:3
    ic=ic+1
    subplot(3,1,i)
    for chan = 8%1:size(itc,3)
    h=imagesc(itc_cond{i}.time,freqs(ff,:),squeeze(itc(i,ff,chan,:,:)))
    hold on
    axis xy
    pause()
    end
    legend(itc_cond{i}.label)
end

end

%%
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
%end


%%
foi = [4 40];
subplot(3,2,[3,4])
%%
p1=semilogx(freqz,squeeze(nanmean(itc_all(:,:),2)),'linewidth',2);
hold on
p2=semilogx(freqz,squeeze(nanmean(itc_all(:,:),2)),'linewidth',2);
set(gca,'fontsize',14,'xtick',[0,4,10,20,30,40,50,60,205])
%xlabel('Freq(Hz)')
ylabel('ITPC')
%ylim([0 0.6])
%xlim([3 300])
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
