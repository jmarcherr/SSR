clear all
assr_startup
ft_defaults

%% Select preprocessed .mat file to process
addpath('./_data/tech_aud')
cd('./_data/tech_aud')
d=dir('*.mat');

%% subject loop
for ss=1:length(d)
    load([d(ss).name]);
    
    cd(rootdir)
    close all
    ff=0;
    fms=[4 207];
    % trial loop
    chan_idx = 0;
    chan=[find(strcmp(data_DG.label,'Cz'))];%
    %% rename data
    data_all= data_DG;fs = data_all.fsample;
    cfg = [];
    cfg.channel = [chan];
    data_all = ft_selectdata(cfg,data_all);
    chan_blink = find(strcmp(data_DG.label,'Fp1'));
    cfg.channel = chan_blink;
    data_blink = ft_selectdata(cfg,data_all);
    
    
    %% divide into trials
    for ff=1:2 % Frequency loop
        fm_freq_oi = fms(ff);
        %init
        epoched_data = [];epoch_var=[];epoch_weighted=[];trial_data = [];
        trial_weights=[];data_tmp=[];
        
        % trigger
        %condvec = [1,3;2,4];
        for kk=1:4 % condition loop
            trial_id = kk;
            
            %init
            epoched_data = [];epoch_var=[];epoch_weighted=[];
            trial_data = [];
            trial_weights=[];data_tmp=[];
            
            trials_oi = length(find(data_all.trialinfo==trial_id));
            
            epoched_data = epoch_data_2(data_all,trial_id,trials_oi);
            epoched_blink = epoch_data_2(data_blink,trial_id,trials_oi);
            
            %% artifact rejection and weighting
            reject_epoch = 0;
            reject_idx = 0;
            
            if fm_freq_oi ==4
                filt_coef = [fm_freq_oi-2 fm_freq_oi+2]; %4 Hz
            else
                filt_coef = [fm_freq_oi-fm_freq_oi/15 fm_freq_oi+fm_freq_oi/15]; % 207 Hz
            end
            
            filt_def = designfilt('bandpassiir','FilterOrder',2, ...
                'HalfPowerFrequency1',filt_coef(1),'HalfPowerFrequency2',filt_coef(2), ...
                'SampleRate',fs,'designmethod','butter');
            
            
            for ii=1:size(epoched_data,1)
                data_tmp = epoched_data(ii,:);
                data_art = epoched_blink(ii,:);
                reject_idx = reject_idx+1;
                if max(abs(data_art)) > 60 % thresholding 80mV
                    reject_epoch = reject_epoch +1;
                    reject_idx = reject_idx-1;
                    %close all
                    %plot(data_art)
                    %pause
                else
                    %filter between filt_coeff
                    filtered_data = filtfilt(filt_def, data_tmp);
                    %compute variance of filtered epoch
                    epoch_var(reject_idx) = var(filtered_data);
                    %Weight the current epoch by its inverse variance (1/var)
                    epoch_weighted(reject_idx,:) = data_tmp./epoch_var(reject_idx);
                end
            end
            
            epoch_var = epoch_var.^(-1); % weights
            
            %Reporting
            clc
            %Number of valid epochs
            valid_epochs = size(epoched_data,1) - reject_epoch;
            fprintf('Number of valid epochs: %d out of %d \n',valid_epochs,...
                size(epoched_data,1))
            %Rejected epochs ratio (in percentage)
            rejt_ratio(ss,kk) = reject_epoch*100/size(epoched_data,1);
            fprintf('Percentage of rejected epochs %.2f %% \n \n', sum(rejt_ratio(ss,:)))
            
            if rejt_ratio > 25
                warning('Epochs rejected ratio over 25% !!!')
            end
            
            
            %%Concatenating into epochs
            
            %number of epochs per trial
            if fm_freq_oi < 207
                epochs_x_trial = 5;
            else
                epochs_x_trial = 5;
            end
            %number of epochs
            idx_tmp = size(epoch_weighted,1);
            %number of valid trials
            nr_trials = floor(idx_tmp/epochs_x_trial);
            %length of each epoch in samples
            sample_x_epoch = length(epoch_weighted);
            
            % generate epochs
            cc=0;
            epoched_data = [];
            
            trial_data = reshape(epoch_weighted(1:nr_trials*epochs_x_trial,:)',...
                sample_x_epoch*epochs_x_trial, nr_trials);
            trial_weights = reshape(epoch_var(1:nr_trials*epochs_x_trial),...
                epochs_x_trial, nr_trials);
            
            
            
            trial_data_cat{kk} = trial_data;
            trial_weights_cat{kk} = trial_weights;
        end
        
        %%
        
        
        
        
        %% Summing over trials
        
        for cond = 1:4 %looping for test-retest
            
            
            D_data = trial_data_cat{cond};
            D_weights = trial_weights_cat{cond};
            
            
            summed_data = sum(D_data,2); % picking random segments
            summed_weights = sum(D_weights,2); % finding weights
            
            idx_sum_weight = 1;
            for idx_sum = 1:sample_x_epoch:length(summed_data)
                data_weighavg(idx_sum:idx_sum + sample_x_epoch-1, 1) = ...
                    summed_data(idx_sum:idx_sum+sample_x_epoch-1)./...
                    summed_weights(idx_sum_weight);
                
                idx_sum_weight = idx_sum_weight + 1;
            end
            M=[];
            M = data_weighavg;
            
            %% %% perform FFT on epoched data
            f_fft = [];
            %FFT
            f_fft = fft(M)/(length(M)/2);
            %Convert to power
            f_fft_pow = abs(f_fft.^2); %
            %Truncate negative freqencies
            f_fft_pow = (f_fft_pow(1:end/2+1));
            %Frequency vector
            f = fs/2*linspace(0,1,length(f_fft_pow));
            
            %% finding peak
            %fc's
            fc_idx = [fm_freq_oi];
            fc = find(f==fc_idx);
            
            bg_freq = [];
            % background (fc +/-3 Hz)
            if fm_freq_oi == 4
                bg_freq = [find(f>=fc_idx-1.5 & f<fc_idx), ...
                    find(f<=fc_idx+1.5 & f>fc_idx)];
            else
                bg_freq = [find(f>=fc_idx-1.5 & f<fc_idx), ...
                    find(f<=fc_idx+1.5 & f>fc_idx)];
            end
            
            nogo_freqs = [4+1/3 4-1/3 207+1/3 207-1/3 207+4 207-4 207+8 207-8];
            for t=1:length(nogo_freqs)
                ng_freq(t)=find(f==nogo_freqs(t));
            end
            bg_freq = setdiff(bg_freq,ng_freq);
            
            %Power at fm_freq_oi
            Psn{ss,ff,cond} = f_fft_pow(fc);
            %Averaged power in surrounding noise bands
            Pn{ss,ff,cond} = mean(f_fft_pow(bg_freq));
            %SNR(dB)
            SNR{ss,ff,cond} = 10*log10(f_fft_pow(fc))-10*log10(mean(f_fft_pow(bg_freq)));
            %F-statistic
            F{ss,ff,cond}=f_fft_pow(fc)/mean(f_fft_pow(bg_freq));
            %Signal power
            Ps{ss,ff,cond} = f_fft_pow(fc)-mean(f_fft_pow(bg_freq));
            F_crit{ff} = finv(0.95,2,2*length(bg_freq));
            
            
            
        end
    end
end

%%

%%
close all
colorcodes = cbrewer('qual','Paired',4)%{'k','g','b','r'};
pow_sig = cellfun(@minus,Psn,Pn,'Un',0);
for ii=1:length(d) %number of subs
    figure(1)
    for f=1:2 %frequency
        subplot(1,2,f)
        for cc=1:4 %condition
            sort = [1,3,2,4];
            if SNR{ii,f,sort(cc)}<db(sqrt(F_crit{ff}-1))
                plot(cc,SNR{ii,f,sort(cc)}, '-o','color',colorcodes(cc,:))
            else
                p(cc)=plot(cc,SNR{ii,f,sort(cc)}, '-o','color',colorcodes(cc,:),'markerfacecolor',colorcodes(cc,:))
            end
            hold on
            %plot(ii,Pn{f,cc,ii},[colorcodes{cc},'x'])
        end
        
        
        plot(1:4,cell2mat(squeeze(SNR(ii,f,sort(:)))),'-','color',[0 0 0 0.1])
        sig = plot(1:4,ones(1,4)*db(sqrt(F_crit{ff}-1)),'--k')
        
        
        box off
        xlim([0.5 4+.5])
        ylim([-15 20])
        set(gca,'xtick',0)
        xlabel('Condition')
        %ylabel('SNR dB')
        
        title([num2str(fms(f)), ' Hz'])
        set(gca,'Fontsize',16)
    end
    
end
hleg = legend([p(1),p(2),p(3),p(4),sig],'25%,40dB', '85%,40dB','25%,75dB','85%,75dB','F_{sig}','location','southeast')
hleg.Box = 'off'





%% level 4Hz 207Hz corr
sublabels = {'1','2','3','4','5','6','7','8','9','10','1retest','3retest'};
figure(2)
%207 difference
%Conditions
%1 = low mod, low level
%2 = low mod, high level
%3 = high mod, low level
%4 = high mod, high level
snr_con = cell2mat(SNR)
p_diff = cell2mat(cellfun(@minus,Ps(:,2,4),Ps(:,2,3),'Un',0)); %207 Hz difference at 75dB
p_4Hz = cell2mat(Ps(:,1,2));
sig_idx = [];
for i=1:length(sublabels)
    if snr_con(i,2,2)<db(sqrt(F_crit{2}-1)) | snr_con(i,2,4)<db(sqrt(F_crit{2}-1)) | snr_con(i,1,2)<db(sqrt(F_crit{2}-1))| snr_con(i,1,4)<db(sqrt(F_crit{2}-1))
        sig_idx(i) = 0;
    else
        sig_idx(i) = 1;
        %textscatter(p_diff(i),p_4Hz(i),{sublabels{i}},'fontsize',14)

    end
    hold on
end

sig_idx([11,12]) = 0;
scatter(p_diff(find(sig_idx)),p_4Hz(find(sig_idx)),'.w')
[r,p]=corr(squeeze(p_diff(find(sig_idx))),squeeze(log10(p_4Hz(find(sig_idx)))))

hold on
textscatter(p_diff(find(sig_idx)),p_4Hz(find(sig_idx)),{sublabels{find(sig_idx)}},'fontsize',14)
lsline
title(['r= ',num2str(r),' p= ',num2str(p)])
xlabel('207Hz(85%-25%) amplitude ')
ylabel('4Hz amplitude')

%plot(p_diff([7,11,12]),p_4Hz([7,11,12]),'o')
%%




