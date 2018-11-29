chan_oi = [4 5 6 9 10    11    39    40    41    44    45    46  38 47];
fid = [4 40];
cmp = flip(cbrewer('div','RdBu',3)); cmp = cmp([1 end],:);
cmp = flip(cbrewer('div','RdBu',2)); cmp = cmp([1 3],:); 
cmp = cbrewer('qual','Pastel1',3); cmp = cmp([1 2],:);
%figure('Position',[100 100 800 400]);
figure('Position',[100 100 600 900]);

p_col = [];
itpc_table = table(SUBJECTS(:),CP(:));
itpc_table.Properties.VariableNames = {'Subject','HI'};
tit = {'Delta - gamma','Gamma'};
titP = {'4 Hz delta - gamma','40 Hz gamma'};

% loop over 4-hz delta-gamma and 40-hz gamma
for cond = 1  :2
    
    
    itpc_subject = [];
    
    % peform the analysis in short time-windows
    for k = 1 : 8
        time_oi = find(time>(0+0.25*(k-1)) & time < ((0.25+0.25*(k-1))));
        xticks{k} = [num2str(0+0.25*(k-1)),' s - ',num2str((0.25+0.25*(k-1))),' s'];
        
        
        itpc_col = [];
        for c = 1 : 2
            itpc_col =[itpc_col squeeze(mean(mean(itpc(chan_oi,time_oi,freq==fid(cond),CP==c-1,cond))))];
        end
        
        % [sbj x hi x time]
        itpc_subject(:,:,k)=itpc_col;
        
        % statistics
        [~,p,~,stats] = ttest2(itpc_col(:,1),itpc_col(:,2));
        p_col(k,cond) = p;

        
        itpct = table(SUBJECTS(:),CP(:),[itpc_col(:,1); itpc_col(:,2)]);
       table_varname = [tit{cond},'_',num2str(k)];
       itpct.Properties.VariableNames = {'Subject','HI',table_varname(table_varname~=' '&table_varname~='-')}

       
           
           
           itpc_table = join(itpc_table,itpct)
        
    end

    

    
    
    
    tmp = mean(itpc_subject,3);
    table_varname = [tit{cond},'_mean'];
    itpct.Properties.VariableNames = {'Subject','HI',table_varname(table_varname~=' '&table_varname~='-')}
    itpc_table = join(itpc_table,itpct)
    
    [~,pall]=ttest2(tmp(:,1),tmp(:,2))
    disp(pall)
    
    subplot(3,2,cond)
    % gamle version
    % errorbar(squeeze(mean(itpc_subject))',squeeze(std(itpc_subject))'/sqrt(40),'LineWidth',2)
    
    for c = 0 : 1
        stdx =  squeeze(std(itpc_subject(:,c+1,:)))./sqrt(40);
        mux = squeeze(mean(itpc_subject(:,c+1,:)));
        if c == 0
            h1=shadedErrorBar(1:8,mux,stdx,'lineprops',{'-','Color',cmp(c+1,:),'MarkerFaceColor',cmp(c+1,:)})
        else
            h2=shadedErrorBar(1:8,mux,stdx,'lineprops',{'-','Color',cmp(c+1,:),'MarkerFaceColor',cmp(c+1,:)})
        end
        hold on
        plot(1:8,mux,'.','Color',cmp(c+1,:),'MarkerSize',12)
    end
    set(gca,'xtick',1:8)
    set(gca,'xticklabel',xticks)
    xtickangle(40)
    
    if cond == 1
        text(1.9,0.54,'*','Fontsize',15)
        text(2.9,0.54,'*','Fontsize',15)
    end
    ylim([0.15 0.6])
    ylabel('ITPC')
    title(titP{cond})
    hleg = legend([h1.mainLine,h2.mainLine],'HI','NH')
    hleg.Box = 'off'
    
    set(gca,'Fontsize',11)
    
end

[p_fdr, p_masked] = fdr( p_col, 0.05);
save('itpc_results','itpc_table')

%%
time = itc_cond{1,1}.time;
time_oi = find(time>0 & time < 2);
%figure('Position',[100 100 800 800]);

%time = itc_cond{1,1}.time;
titA ={'HI','NH'};
titB ={'delta-gamma','gamma'};
hi=0;
for k = 1 :2
        subplot(3,2,k+(hi)*2+2)
        imagesc(time(time_oi),freq,squeeze(mean(mean(itc(:,time_oi,:,CP==hi,k),1),4))',[0.05 0.6]); axis xy;
        cmp = flip(cbrewer('div','RdBu',100));
        colormap(cmp);
        ylabel('Frequency (Hz)')
        xlabel('Time (s)')
        title({[titA{hi+1},' ' titB{k}]})
        ylim([2 50])
        if k==1
            yticks = [4 10 20 30 40 50];
        else
            yticks = [10 20 30 40 50];
        end
        set(gca,'ytick',[yticks])
        set(gca,'yticklabel',([yticks]))
        
        set(gca,'Fontsize',11)
end



