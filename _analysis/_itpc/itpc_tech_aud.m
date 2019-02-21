clear all
% ITPC analysis
assr_startup    % add your paths to data and scripts
ft_defaults     % fieldtrip defaults
%% Select preprocessed .mat file to process
addpath('./_data/tech_aud/itpc')
cd('./_data/tech_aud/itpc')
d=dir('*.mat');

%% subject loop
for ss=1:length(d)
    load([d(ss).name]);
    data_all = data_DG;
    cd(rootdir)
    close all
    ff=0;
    fms=[4 207];
    % trial loop
    chan_idx = 0;
    chan=[find(strcmp(data_DG.label,'Cz'))];%
%%%%%%%%%%%
%% Threshold denoising

rejectidx = [];
trigs=[1,2,3,4];
for tt=1:length(trigs) % loop over triggers
    idx = find(data_all.trialinfo==trigs(tt));
    for ii=1:length(idx)
        if max(abs(mean(data_all.trial{idx(ii)}(1:2,:))))>120 % 120 mv at frontal channels (Fp1,Fp2)
            rejectidx = [rejectidx idx(ii)];
        else
        end
    end
    
end
disp([num2str((length(rejectidx)/length(data_all.trialinfo))*100),'% of trials rejected'])
%% Phase coherence (single frequency)

foihz =[4;207]; % Frequencies of interest
for ff  = 1 : size(foihz,1)
    for k=1:length(trigs)
        %% ITPC
        idx=setdiff(find(data_all.trialinfo==trigs(k)),rejectidx)
        idx_tmp=idx;
        cfg = [];
        cfg.method = 'wavelet';
        cfg.toi    = -3:0.01:5;
        cfg.foi    = foihz(ff,:);
        cfg.trials       = idx_tmp;
        cfg.output = 'fourier';
        cfg.channel      = {'all', '-Status'};
        cfg.pad = 'nextpow2';
        if ff>1
            cfg.width = 500;
        else
            cfg.width =12;%<- like dtu /hvidovre paper
        end
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

%%
% ITC
%Gather all
itc_all=[];
itc = [];
freqs = [];
for ff=1:2
    for kk=1:length(trigs)
        itc(kk,ff,:,:,:) = itc_cond{kk,ff}.itpc; %cond x freq x chan x time
        freqs(ff,:)= itc_cond{kk,ff}.freq;       % frequency vector
    end
end
itpc_sub{ss} = itc;
end
%%
close all
toi = find(itc_cond{1}.time>=1 & itc_cond{1}.time<1.5)
colorcodes = cbrewer('qual','Paired',4);
sort = [1,3,2,4];
fids = [{'4Hz'},{'207Hz'}];
for ff=1:2
    figure(ff)
for kk=1:4
for ss=1:12

itpc_all(ss,:,:,:,:) = cell2mat(itpc_sub(ss));
plot(itc_cond{1}.time,squeeze(mean(itpc_all(ss,sort(kk),ff,chan,:),1)),'color',[colorcodes(kk,:) 0.5])
hold on

end
p(kk)=plot(itc_cond{1}.time,squeeze(mean(itpc_all(:,sort(kk),ff,chan,:),1)),'color',colorcodes(kk,:),'linewidth',4)
itpc_avg_sub{ff,kk} = squeeze(mean(itpc_all(:,kk,ff,chan,toi),5))
end

xlim([-0.5 3.5])
ylim([0 0.6])
plot(zeros(1,10),linspace(0,1,10),'--','color',[0.5 0.5 0.5 0.5])
plot(ones(1,10)*3,linspace(0,1,10),'--','color',[0.5 0.5 0.5 0.5])
hleg = legend([p(1),p(2),p(3),p(4)],'25%,40dB', '85%,40dB','25%,75dB','85%,75dB','F_{sig}','location','northeast')
hleg.Box = 'off'

ylabel('ITPC')
xlabel('Time (s)')

set(gca,'Fontsize',16)
title(fids{ff})
box off
set(gcf,'position',[421 646 388 420])

end



%%
%% 2d plot
close all
tit = {'Burst AM','Double SAM','Continuous double SAM'};
choi = find(strcmp(data_all.label,'Cz')); % Potential between P9+10 and Cz
for ff=1:2
    for i=1:4
        subplot(3,1,i)
        h=plot(itc_cond{i}.time,squeeze(itc(i,ff,choi,:)));
        title(cell2mat(tit(i)));
        hold on
        if i==3
            xlabel('Time (s)');
        end
        ylabel('ITPC');
    end
    xlim([-1 4])
    ylim([0 0.4])
end


legend('4Hz ITPC','207Hz ITPC','location','best')
set(gcf,'Position',[467 382 310 535])