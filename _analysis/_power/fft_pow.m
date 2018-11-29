clear all
assr_startup
ft_defaults

%% Select preprocessed .mat file to process
cd(datadir)
load(['subj-ando90_exp-deaff_mod_depth_4m_pilot_ear-Right_lvl-81dB_mod-85.mat']);
cd(rootdir)
%% rename data
data_all= data_DG;


%% Select trials based on triggers (see preprocessing)
idx=find(data_all.trialinfo==1)
idx_tmp=idx;%(1:300);
% how many epochs to group together
epoch = 300;
% generate epochs
cc=0;
epoched_data = [];
for ii = 1:epoch:length(idx_tmp)
    cc=cc+1;
    cfg = [];
    cfg.trials = idx_tmp(ii:ii+epoch-1);
    data_tmp = ft_selectdata(cfg,data_all);
    epoched_data(cc,:,:) = cell2mat(data_tmp.trial(:)'); % epoch x chan x time

end
 

    
%% perform FFT on epoched data
f_fft = [];fs = data_all.fsample;
M = [];
for kk=1:size(epoched_data,1) % run through epochs
    M=squeeze(epoched_data(kk,:,:));
    f_fft(kk,:,:) = fft(M,[],2)./(size(M,2)./2); % length of vector
    %f_fft(kk,:,:) = fft(M,2.^nextpow2(length(M)),2); %nextpow2
end
f_fft_mean = squeeze(nanmean(f_fft,1));
f_fft_pow = abs(f_fft_mean.^2);
f = linspace(0,fs-1/fs,length(f_fft_pow)); % frequency vector


%%
figure(1)
plot(f,db(f_fft_pow(find(strcmp(data_tmp.label,'Cz')),:)));
xlim([0 110])



