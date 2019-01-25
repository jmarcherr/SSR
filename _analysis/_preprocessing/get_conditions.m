%Get conditions

clear all
load(['subject_info.mat'])
%%
Data.conditions
idx=0;
for ii=1:length(Data)
    if ~isempty(Data(ii).conditions)
        idx = idx+1;
        cond{idx} = Data(ii).conditions;
    end
end

clear c
for cc=1:length(cond)
    tmp = cond{cc}
    for tt=1:4
        if tmp(tt,1) ==1 % low level
            if tmp(tt,2) ==1 %low mod
                c(cc,tt) = 1;
            else
                c(cc,tt) = 2; %low level high mod
            end
        elseif tmp(tt,1) ==2 %high level
            if tmp(tt,2) ==1 %low mod
                c(cc,tt) =3; %high level low mod
            else
                c(cc,tt) =4; %high level high mod
            end
        end
    end
    
end


%1=low mod, low level
%2=low mod, high level
%3=high mod, low level
%4=high mod, high level



%%
