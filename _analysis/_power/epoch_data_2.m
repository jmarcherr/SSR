function [data_out] = epoch_data_2(data_in,trial_id,nr_epochs)
%Epoch data

%More info following

if ~isstruct(data_in)
    error('Data not in right format');
else
    clc
    disp('Epoching data');
end
if length(trial_id)>1
    idx = find(data_in.trialinfo==trial_id(1) | data_in.trialinfo==trial_id(2))
else
    idx=find(data_in.trialinfo==trial_id)
end
idx_tmp=idx(nr_epochs);

% generate epochs
epoched_data = [];
for ii = 1:length(idx_tmp)
    cfg = [];
    cfg.trials = idx_tmp(ii);
    %data_tmp = ft_selectdata(cfg,data_in);
    epoched_data(ii,:,:) = cell2mat(data_in.trial(idx_tmp(ii))); % epoch x chan x time
    clc
    disp(['trial ',num2str(ii), ' of ', num2str(length(idx_tmp)), ' done.'])
end

if size(epoched_data,1)>1 %mean over channels
    epoched_data = mean(epoched_data,2);
end
data_out = squeeze(epoched_data);

end
