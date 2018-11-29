function [data_clean,t_count,art_count] = thr_rjct(data_erp,thres,pad,whole_trial,plot)

fs = data_erp.fsample;
art=[];
data_clean = data_erp;
t_count = 0;
art_count = 0;

for tt=1:length(data_erp.trial)
    for cc = 1:size(data_erp.trial{1},1)
        art{tt,cc} = find(abs(data_erp.trial{tt}(cc,:))>thres);
        
        if ~isempty(art{tt,cc})
            if whole_trial
                art_time = 1:length(data_erp.trial{tt});
            else
            % Padding for artefacts
            if min(art{tt,cc})<=round(pad*fs)+1 % less than pad sec
                art_start = 1;%min(art{tt,cc});
            else
                art_start = min(art{tt,cc})-(round(pad*fs)-1);
            end
            if max(art{tt,cc})>=(length(data_erp.trial{tt}(cc,:))-round(pad*fs))
                art_end = length(data_erp.trial{tt}(cc,:));
            else
                art_end = max(art{tt,cc})+round(pad*fs);
            end
            
            % Padded artefact time 
            art_time = art_start:art_end;
            end
            if plot && cc==38
                plot(data_erp.time{tt}',data_erp.trial{tt}(cc,:))
                hold on
                plot(data_erp.time{tt}(art_time)',data_erp.trial{tt}(cc,art_time),'r')
                pause
                close all
            end
            
            

            
            % 
            if art_time(end)>129
             %   disp('error')
            end
            data_clean.trial{tt}(:,art_time) = 0;
            %sec_reject = [sec_reject + length(art_time)/fs]
        end
        
    end
    
    
    % counting artefacts
    if (sum(cat(2,art{tt,:},1))-1)
        t_count = t_count+1;
        art_count = art_count +length(art_time);
    end
end

% output the results
disp(['found ' num2str(art_count) ' artefacts:' ' filled ' num2str(t_count) ' trials with nans'])
