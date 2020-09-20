function meanChoiceSelectivityAcrossSessions

% this functions calculate mean choice selectivity across sessions per
% animal, takes about ~5 minutes to finish. 321.426

%%  Choice selectivity analysis
dffFilePath = 'D:\JenHau\siniscalchi2019\Learning\analysis';
fpath='D:\JenHau\siniscalchi2019\Learning\longitudinal\figures\Choice selectivity';
animalList = [{'M52'};{'M53'};{'M54'};{'M55'};{'M56'}];

% Calculate PSTH values first!
tic
params=[];
params.xtitle = 'Time from stimulus (s)';
params.window = [-2:0.25:6.5];
params.minNumTrial = 5;

for animalID = 1: numel(animalList)
    curr_animal = animalList{animalID};
    temp_an = dir(fullfile(dffFilePath,['*',curr_animal,'*']));
    fields = [ {'choicesel_hit'},{'choicesel_err'},{'choicesel_left'},{'choicesel_right'}]; 
    % hit = L vs. R; err = L_err vs. R_err; left = L hit vs. err; right = R hit vs. err 
    
    for ses=1: size(temp_an,1)
        load(fullfile(temp_an(ses).folder,temp_an(ses).name,'dff.mat'));
        load(fullfile(temp_an(ses).folder,temp_an(ses).name,'beh.mat'));
        params.trigTime = trialData.cueTimes;
        
        % Calculate PSTH for all cells
        psth_left = [];
        psth_right = [];
        psth_left_err = [];
        psth_right_err = [];       
        choicesel_hit = [];
        choicesel_err = [];
        choicesel_left = [];
        choicesel_right = [];
        
        for j=1:numel(cells.dFF) %notes: same as meanPSTHAcrossSessions.m
                fieldname={'sound','upsweep','left','hit'}; trialMask = getMask(trials,fieldname);
                psth_left{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
                fieldname={'sound','downsweep','right','hit'}; trialMask = getMask(trials,fieldname);
                psth_right{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);               
                fieldname={'sound','downsweep','left','err'}; trialMask = getMask(trials,fieldname);
                psth_left_err{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
                fieldname={'sound','upsweep','right','err'}; trialMask = getMask(trials,fieldname);
                psth_right_err{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
                
                %calculate choice selectivity = (A-B)/(A+B)
                choicesel_hit{j} = calc_selectivity(psth_left{j},psth_right{j});
                choicesel_err{j} = calc_selectivity(psth_left_err{j},psth_right_err{j});
                choicesel_left{j} = calc_selectivity(psth_left{j},psth_left_err{j});
                choicesel_right{j} = calc_selectivity(psth_right{j},psth_right_err{j});  
        end
        % keep the psth values for each session
        dat(animalID).([fields{1},num2str(ses)]) = choicesel_hit;
        dat(animalID).([fields{2},num2str(ses)]) = choicesel_err;
        dat(animalID).([fields{3},num2str(ses)]) = choicesel_left;
        dat(animalID).([fields{4},num2str(ses)]) = choicesel_right;

        % plot average choice selectivity for each session
        nCell = numel(cells.dFF);
        h = figure('Name',(curr_animal),'NumberTitle','off');
        subplot(2,1,1); hold on
        val_sig = nan(nCell, numel(choicesel_hit{1}.signal));
        for j=1:nCell
            val_sig(j,:) =choicesel_hit{j}.signal;
        end
        t = choicesel_hit{1}.t;
        plot(t,nanmean(val_sig),'color','k');
        title ('Hit Trials')
        box off
        xlabel ('Time (seconds)')
        ylabel ('Choice selectivity')
        
        subplot(2,1,2); hold on
        val_sig = nan(nCell, numel(choicesel_hit{1}.signal));
        for j=1:nCell
                val_sig(j,:) =choicesel_err{j}.signal;
        end
        t = choicesel_hit{1}.t;
        plot(t,nanmean(val_sig),'color','k');
        title ('Error Trials')
        box off
        xlabel ('Time (seconds)')
        ylabel ('Choice selectivity')
        
        saveas(figure(h),fullfile(fpath,[(curr_animal),'session',num2str(ses)]),'fig')
        saveas(figure(h),fullfile(fpath,[(curr_animal),'session',num2str(ses)]),'png')
        close all
    end
end
toc


%% plot choice selectivity 
fpath='D:\JenHau\siniscalchi2019\Learning\longitudinal\figures\Choice selectivity';
params.xtitle = 'Time from stimulus (s)';
% params.sortParam=[0 6.5];  %sort the cells based on amplitude of selectivity in this time period
params.colorRange=[-0.8 0.8];
tlabel = ['Choice selectivity'];

% hit vs. err
for animalID = 1: numel(animalList)
    curr_animal = animalList{animalID};
    h = figure('Name',(curr_animal),'NumberTitle','off');
    params.xtitle = 'Time from stimulus (s)';    
    for ses =1:5
        subplot(2,5,ses)       
        temp = dat(animalID).(['choicesel_left',num2str(ses)]);
        nCells = numel(temp);
        params.sortParam=[1:nCells]; %do not sort the cells, list from #1 to #nCell
        plot_selectivity(temp,params.sortParam,tlabel,params.xtitle,params.colorRange);
        title(['Session: ',num2str(ses, '%.2d')]);
        
        subplot(2,5,ses+5)
        temp = dat(animalID).(['choicesel_right',num2str(ses)]);
        plot_selectivity(temp,params.sortParam,tlabel,params.xtitle,params.colorRange);
        title(['Session: ',num2str(ses, '%.2d')]);
    end
    sgtitle('Choice selectivity, left vs. right');
    
%     %make a color scale bar
%     subplot(3,20,60);
%     image(0,linspace(colorRange(1),colorRange(2),100),linspace(colorRange(1),colorRange(2),100)','CDataMapping','scaled');
%     colormap(colors);
%     caxis([colorRange(1) colorRange(2)]);
%     title(['(A-B)/(A+B)']);
%     set(gca,'YDir','normal');
%     set(gca,'XTick',[]);
    
    saveas(figure(h),fullfile(fpath,(curr_animal)),'fig')
    saveas(figure(h),fullfile(fpath,(curr_animal)),'png')
end

%% %%
% colors=cbrewer('div','RdBu',256);
% colors=flipud(colors);
% colorRange=[-0.8 0.8];
% 
% t=temp{1}.t;
% nCells=numel(temp);  
% 
% pref=[];
% for j=1:nCells
%     pref(:,j)=temp{j}.signal;
% end
% 
%   figure;
% 
% % if numel(sortParam) == 2
% 
%     tIdx=[max([sum(t<=sortParam(1)) 1]):sum(t<=sortParam(2))];  %index should start from at least value of 1
% 
%     negPrefCells = find(nanmean(pref(tIdx,:),1)<0);  %determine sign of preference
%     posPrefCells = find(nanmean(pref(tIdx,:),1)>=0);
%     
%     for j=1:numel(negPrefCells)     % sort by center of mass (mass should be all positive)
%         mass = pref(tIdx,negPrefCells(j));
%         mass(mass>0) = 0;
%         com_neg(j) = -sum(t(tIdx).*mass)/sum(mass);
%     end
%     [~,neg_idxSort]=sort(com_neg);
%     for j=1:numel(posPrefCells)
%         mass = pref(tIdx,posPrefCells(j));
%         mass(mass<0) = 0;
%         com_pos(j) = sum(t(tIdx).*mass)/sum(mass);
%     end
%     [~,pos_idxSort]=sort(com_pos);
%     cellOrder = [negPrefCells(neg_idxSort) posPrefCells(pos_idxSort)];
%     
% % elseif numel(sortParam) == nCells
%     %sort by a specified order
% %     cellOrder = sortParam; 
% % else    
% %     error('Error with the sortParam input for the plot_selectivity() function.');
% % end
% 
% %plot in pseudocolor
% % subplot(1,3,1);
% image(t,1:nCells,pref(:,cellOrder)','CDataMapping','scaled');
% hold on; plot([0 0],[0 nCells+1],'w');
% colormap(colors);
% caxis([colorRange(1) colorRange(2)]);      %normalize dF/F heatmap to max of all conditions
% ylabel('Cells');
% % xlabel(xtitle);
% title({tlabel;['A=' temp{1}.input1_label];['B=' temp{1}.input2_label]});
% 

