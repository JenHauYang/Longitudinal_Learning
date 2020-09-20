% calculates how many choice/outcome selective cells, takes ~ an hour

%% CHOICE AND OUTCOME: Multiple linear regression  - choice and reward and their interaction
dffFilePath = 'D:\JenHau\siniscalchi2019\Learning\analysis';
animalList = [{'M52'};{'M53'};{'M54'};{'M55'};{'M56'}];

tic
for animalID = 1: numel(animalList)
    curr_animal = animalList{animalID};
    temp_an = dir(fullfile(dffFilePath,['*',curr_animal,'*']));
    fields = [ {'choicesel_hit'},{'choicesel_err'},{'choicesel_left'},{'choicesel_right'}];
    sel_all = nan(1,size(temp_an,1));
    
    for ses=1: size(temp_an,1)
        load(fullfile(temp_an(ses).folder,temp_an(ses).name,'dff.mat'));
        load(fullfile(temp_an(ses).folder,temp_an(ses).name,'beh.mat'));
        
        params=[];
        params.trigTime = trialData.cueTimes;
        % first predictor is choice; dummy-code: left=-1, right=1, miss=NaN
        params.choiceEvent=NaN(size(trials.left));
        params.choiceEvent(trials.left) = -1;
        params.choiceEvent(trials.right) = 1;
        % second predictor is outcome; dummy-code: reward=1, miss=0
        params.outcomeEvent=NaN(size(trials.hit));
        params.outcomeEvent(trials.hit) = 1;
        params.outcomeEvent(trials.err) = 0;
        
        params.xtitle = 'Time from stimulus (s)';
        params.window = [-2:0.5:6.5];
        params.nback = 2;       %how many trials back to regress against
        params.interaction = true; %consider interaction terms (our data do not have enough trials)
        params.pvalThresh = 0.01;   %p-value for coefficient be considered significant
        
        % only perform analysis on trials with a response (when animal is engaged)
        fieldname={'left','right'};
        % fieldname={'hit','err'};
        trialMask = getAnyMask(trials,fieldname);
        
        for j=1:numel(cells.dFF)
            reg_cr{j}=linear_regr( cells.dFF{j}, cells.t, [params.choiceEvent params.outcomeEvent], params.trigTime, trialMask, params );
        end
        % tlabel={'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'};
 
        % which cell was choice selective
        timeIdx=sum(0>reg_cr{1}.regr_time);   %find index associated with time = 0 s
        
        nCells = numel(reg_cr);
        choice_mod=false(nCells,1);
        choice1_mod=false(nCells,1);
        choice2_mod=false(nCells,1);
        outcome_mod=false(nCells,1);
        outcome1_mod=false(nCells,1);
        outcome2_mod=false(nCells,1);
        interaction_mod=false(nCells,1);
        interaction1_mod=false(nCells,1);
        interaction2_mod=false(nCells,1);
        
        for j=1:nCells
            if sum(reg_cr{j}.pval(timeIdx:end,2)<params.pvalThresh) >= 5 %2 for C(n)
                choice_mod(j)=true;
            end
            if sum(reg_cr{j}.pval(timeIdx:end,3)<params.pvalThresh) >= 5 %3 for C(n-1)
                choice1_mod(j)=true;
            end
            if sum(reg_cr{j}.pval(timeIdx:end,4)<params.pvalThresh) >= 5 %4 for C(n-2)
                choice2_mod(j)=true;
            end

            if sum(reg_cr{j}.pval(timeIdx:end,5)<params.pvalThresh) >= 5 %5 for R(n)
                outcome_mod(j)=true;
            end
            if sum(reg_cr{j}.pval(timeIdx:end,6)<params.pvalThresh) >= 5 %6 for R(n-1)
                outcome1_mod(j)=true;
            end
            if sum(reg_cr{j}.pval(timeIdx:end,7)<params.pvalThresh) >= 5 %7 for R(n-2)
                outcome2_mod(j)=true;
            end
            
            if sum(reg_cr{j}.pval(timeIdx:end,8)<params.pvalThresh) >= 5 %8 for C(n)xR(n)
                interaction_mod(j)=true;
            end
            if sum(reg_cr{j}.pval(timeIdx:end,9)<params.pvalThresh) >= 5 %9 for C(n-1)xR(n-1)
                interaction1_mod(j)=true;
            end
            if sum(reg_cr{j}.pval(timeIdx:end,10)<params.pvalThresh) >= 5 %10 for C(n-2)xR(n-2)
                interaction2_mod(j)=true;
            end
        end        
        mod = choice_mod | interaction_mod;
        mod1 = choice1_mod | interaction1_mod;
        mod2 = choice2_mod | interaction2_mod;
        sel_choice(ses) = sum(mod);
        sel_choice1(ses) = sum(mod1);
        sel_choice2(ses) = sum(mod2);
        sel_meanchoice(ses) = mean(mod);
        sel_meanchoice1(ses) = mean(mod1);
        sel_meanchoice2(ses) = mean(mod2);
        
        mod = outcome_mod | interaction_mod;
        mod1 = outcome1_mod | interaction1_mod;
        mod2 = outcome2_mod | interaction2_mod;
        sel_outcome(ses) = sum(mod);
        sel_outcome1(ses) = sum(mod1);
        sel_outcome2(ses) = sum(mod2);
        sel_meanoutcome(ses) = mean(mod);
        sel_meanoutcome1(ses) = mean(mod1);
        sel_meanoutcome2(ses) = mean(mod2);
        
    end    
    selectivity(animalID).choice = sel_choice;
    selectivity(animalID).choice1 = sel_choice1;
    selectivity(animalID).choice2 = sel_choice2;
    selectivity(animalID).meanchoice = sel_meanchoice;
    selectivity(animalID).meanchoice1 = sel_meanchoice1;
    selectivity(animalID).meanchoice2 = sel_meanchoice2;
    selectivity(animalID).outcome = sel_outcome;
    selectivity(animalID).outcome1 = sel_outcome1;
    selectivity(animalID).outcome2 = sel_outcome2;
    selectivity(animalID).meanoutcome = sel_meanoutcome;
    selectivity(animalID).meanoutcome1 = sel_meanoutcome1;
    selectivity(animalID).meanoutcome2 = sel_meanoutcome2;
end
toc    
    
% plot numbers of choice selective cells across sessions
fpath='D:\JenHau\siniscalchi2019\Learning\longitudinal\figures';
analysisFilePath = 'D:\JenHau\siniscalchi2019\Learning\analysis';
setup_figprop;

mice_sel = struct2cell(selectivity);
choice_sel = cell2mat(squeeze(mice_sel(1,:,:)));
choice1_sel = cell2mat(squeeze(mice_sel(2,:,:)));
choice2_sel = cell2mat(squeeze(mice_sel(3,:,:)));
meanchoice_sel = cell2mat(squeeze(mice_sel(4,:,:)));
meanchoice1_sel = cell2mat(squeeze(mice_sel(5,:,:)));
meanchoice2_sel = cell2mat(squeeze(mice_sel(6,:,:)));
outcome_sel = cell2mat(squeeze(mice_sel(7,:,:)));
outcome1_sel = cell2mat(squeeze(mice_sel(8,:,:)));
outcome2_sel = cell2mat(squeeze(mice_sel(9,:,:)));
meanoutcome_sel = cell2mat(squeeze(mice_sel(10,:,:)));
meanoutcome1_sel = cell2mat(squeeze(mice_sel(11,:,:)));
meanoutcome2_sel = cell2mat(squeeze(mice_sel(12,:,:)));
beh = [1:size(temp_an,1)];

h = figure;
for k=1:6
    subplot(2,3,k);hold on
    if k ==1
        temp_fig = choice_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'C(n) selective cells';
        yrange = [0 100];
    elseif k ==2
        temp_fig = choice1_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'C(n-1) selective cells';
        yrange = [0 100];
    elseif k ==3
        temp_fig = choice2_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'C(n-2) selective cells';
        yrange = [0 100];
    elseif k ==4
        temp_fig = outcome_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'R(n) selective cells';
        yrange = [0 150];
    elseif k ==5
        temp_fig = outcome1_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'R(n-1) selective cells';
        yrange = [0 150];
    elseif k ==6
        temp_fig = outcome2_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'R(n-2) selective cells';
        yrange = [0 150];  
    end    
        plot(beh,temp_fig, '.:','color',[0.5 0.5 0.5]); hold on
        plot(beh, mean_temp_fig,'s-','color','k')
        xlim([0 size(temp_an,1)+1]); box off
        xticks([1:5])
        ylim(yrange)
        xlabel('Session')
        ylabel(tName)
end
saveas(figure(h),fullfile(fpath,['selective_cells.fig']),'fig')
saveas(figure(h),fullfile(fpath,['selective_cells.png']),'png')

h = figure;
for k=1:6
    subplot(2,3,k);hold on
    if k ==1
        temp_fig = 100*meanchoice_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'C(n) selective cells (%)';
        yrange = [0 40];
    elseif k ==2
        temp_fig = 100*meanchoice1_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'C(n-1) selective cells (%)';
        yrange = [0 40];
    elseif k ==3
        temp_fig = 100*meanchoice2_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'C(n-2) selective cells (%)';
        yrange = [0 40];
    elseif k ==4
        temp_fig = 100*meanoutcome_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'R(n) selective cells (%)';
        yrange = [0 40];
    elseif k ==5
        temp_fig = 100*meanoutcome1_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'R(n-1) selective cells (%)';
        yrange = [0 40];
    elseif k ==6
        temp_fig = 100*meanoutcome2_sel';
        mean_temp_fig = nanmean(temp_fig,2);
        tName= 'R(n-2) selective cells (%)';
        yrange = [0 40];  
    end    
        plot(beh,temp_fig, '.:','color',[0.5 0.5 0.5]); hold on
        plot(beh, mean_temp_fig,'s-','color','k')
        xlim([0 size(temp_an,1)+1]); box off
        xticks([1:5])
        ylim(yrange)
        xlabel('Session')
        ylabel(tName)
end
saveas(figure(h),fullfile(fpath,['selective_cells_percentage.fig']),'fig')
saveas(figure(h),fullfile(fpath,['selective_cells_percentage.png']),'png')
       
%% %% Plot cue-aligned dF/F for a representative cell
%     params=[];
%     params.trigTime = trialData.cueTimes;
%     params.xtitle = 'Time from sound cue (s)';
%     params.window = [-3:0.5:7];
%     params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
%     params.CI = 0.9;  %confidence interval
%     params.minNumTrial = 0; %only calc PSTH if there is this number of trials
% 
% % plot different choices in same panels
%     psth_panel=[];
%     i=1;j=6; % a representative cell
%     k=1;
%     fieldname=[];
%     fieldname{1}={'sound','upsweep','left','hit'}; col{1}='r'; linstyle{1}='-';
%     fieldname{2}={'sound','downsweep','right','hit'}; col{2}='b'; linstyle{2}='-';
%         for kk=1:numel(fieldname)
%             trialMask = getMask(trials,fieldname{kk});
%             psth_panel(k).sig{kk} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
%             psth_panel(k).col{kk} = col{kk};
%             psth_panel(k).linstyle{kk} = linstyle{kk};
%         end
%     tlabel = ['Cell ' int2str(j)];
%     plot_psth(psth_panel,tlabel,params.xtitle);
%     % print(gcf,'-dpng',['cell' int2str(j) '-choice']);
%     % saveas(gcf, ['cell' int2str(j) '-choice'], 'fig');
%
% Plot fluorescence data, for discrim tasks with no flexibility
% process data files
% %i=2; j=38; %cell example 1
% i=5; j=38; %cell example 2
% %i=1; j=9;  %cell example 3
% 
% disp(['Processing file ' int2str(i) ' out of ' int2str(numel(expData)) '.']);
% 
% % setup/create subdirectories to save analysis and figures
% if isfield(expData(i),'onefolder')   %all the log files saved in one folder
%     temp = sprintf('%04d',i);
%     savematpath = fullfile(dirs.analysis,expData(i).sub_dir,temp);
%     savefluofigpath = fullfile(dirs.analysis,expData(i).sub_dir,temp,'figs-fluo');
% else
%     savematpath = fullfile(dirs.analysis,expData(i).sub_dir);
%     savefluofigpath = fullfile(dirs.analysis,expData(i).sub_dir,'figs-fluo');
% end
% if ~exist(savefluofigpath,'dir')
%     mkdir(savefluofigpath);
% end
% 
% % load the saved behavioral analysis (from start_beh.m)
% cd(savematpath);
% load('beh.mat');
% load('dff.mat');
% 
% cd(savefluofigpath);
% Plot different outcomes in same panels
% 
% psth_panel=[];
% for k=1:3
%     fieldname=[];
%     if k==1 %panel 1
%         fieldname{1}={'sound','upsweep','left','hit'}; col{1}='r'; linstyle{1}='-';
%         fieldname{2}={'sound','upsweep','left','doublereward'}; col{2}='r'; linstyle{2}=':';
%         fieldname{3}={'sound','upsweep','left','omitreward'}; col{3}='r'; linstyle{3}='--';
%     elseif k==2 %panel 2
%         fieldname{1}={'sound','downsweep','right','hit'}; col{1}='b'; linstyle{1}='-';
%         fieldname{2}={'sound','downsweep','right','doublereward'}; col{2}='b'; linstyle{2}=':';
%         fieldname{3}={'sound','downsweep','right','omitreward'}; col{3}='b'; linstyle{3}='--';
%     elseif k==3 %panel 3
%         fieldname{1}={'sound','hit'}; col{1}='k'; linstyle{1}='-';
%         fieldname{2}={'sound','doublereward'}; col{2}='k'; linstyle{2}=':';
%         fieldname{3}={'sound','omitreward'}; col{3}='k'; linstyle{3}='--';
%     end
%     for kk=1:numel(fieldname)
%         trialMask = getMask(trials,fieldname{kk});
%         psth_panel(k).sig{kk} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
%         psth_panel(k).col{kk} = col{kk};
%         psth_panel(k).linstyle{kk} = linstyle{kk};
%     end
% end
% tlabel = ['Cell ' int2str(j)];
% 
% plot_psth(psth_panel,tlabel,params.xtitle);
% print(gcf,'-dpng',['cell' int2str(j) '-outcome']);
% saveas(gcf, ['cell' int2str(j) '-outcome'], 'fig');
