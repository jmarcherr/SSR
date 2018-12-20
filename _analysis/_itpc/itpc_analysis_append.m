clear all
% ITPC analysis
assr_startup    % add your paths to data and scripts
ft_defaults     % fieldtrip defaults

load(['pilot3_2xsamtone+oldx10min_itpc.mat']); % preprocessed data
data_all= data_DG;

%%%%%%%%%%%
%% Threshold denoising

rejectidx = [];
trigs=[1,2,3];
for tt=1:3 % loop over triggers
    idx = find(data_all.trialinfo==trigs(tt));
    for ii=1:length(idx)
        if max(abs(mean(data_all.trial{idx(ii)}(1:2,:))))>120 % 120 mv at frontal channels (Fp1,Fp2)
            rejectidx = [rejectidx idx(ii)];
        else
        end
    end 

end
    disp([num2str((length(rejectidx)/length(data_all.trialinfo))*100),'% of trials rejected'])
%% Phase coherence

foihz =[200:0.1:210]% [0:10];%;[1:0.1:10]
for ff  = 1 : size(foihz,1)   
    for k=1:length(trigs)
        %% ITPC
        idx=setdiff(find(data_all.trialinfo==trigs(k)),find(rejectidx==1))
        idx_tmp=idx;
        cfg = [];
        cfg.method = 'wavelet';
        cfg.toi    = -3:0.01:5;
        cfg.foi    = foihz(ff,:);
        cfg.trials       = idx_tmp;
        cfg.output = 'fourier';
        cfg.channel      = {'all', '-Status'};
        cfg.pad = 'nextpow2'%maxperlen';
        cfg.width = 500% %; %12<- like dtu /hvidovre paper
        %cfg.gwidth = 10;
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
freqs = [];
for ff=1:1
    for kk=1:3
        itc(kk,ff,:,:,:) = itc_cond{kk,ff}.itpc;
        freqs(ff,:)= itc_cond{kk,ff}.freq;
    end
end

close all

figure(2)

tit = {'Burst AM','Double SAM','Continuous double SAM'};

ic=0;
cmp = cbrewer('seq','Blues',100,'cubic');
for ff=1%:2
for i=1:3
    ic=ic+1
    subplot(3,1,i)
    for chan = 8 %Cz
    h=imagesc(itc_cond{i}.time,freqs(ff,:),squeeze(itc(i,ff,chan,:,:)),[0.1 0.25]);
    title(cell2mat(tit(i)));
    hold on
    axis xy
    colormap(cmp)
    c=colorbar;
    c.Label.String = 'ITPC';
    if i==3
        xlabel('Time (s)');
    end
        ylabel('Frequency');
    end
    xlim([-1 4])
    %ylim
end

end

set(gcf,'Position',[467 382 310 535])

%% 2d plot
ic=0;
cmp = cbrewer('seq','Blues',100,'cubic');
fid = find(freqs==207);
for ff=1%:2
for i=1:3
    ic=ic+1
    subplot(3,1,i)
    for chan = 8%1:size(itc,3)
    h=plot(itc_cond{i}.time,squeeze(itc(i,ff,chan,fid,:)));
    title(cell2mat(tit(i)));
    hold on
    axis xy
    colormap(cmp)
    %c=colorbar;
    %c.Label.String = 'ITPC';
    if i==3
        xlabel('Time (s)');
    end
        ylabel('ITPC');
    end
    xlim([-1 4])
    ylim([0 0.4])
end

end

set(gcf,'Position',[467 382 310 535])


