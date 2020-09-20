% notes:

% ccmap =['k','b','r','g','y']; 

%% from choice_selectivity_across_session.m
%
% for animalID = 1: numel(animalList)
%     curr_animal = animalList{animalID};
%     temp_an = dir(fullfile(dffFilePath,['*',curr_animal,'*']));    
%     for ses=1: size(temp_an,1)
%         % load data
%         load(fullfile(temp_an(ses).folder,temp_an(ses).name,'dff.mat'));
%         load(fullfile(temp_an(ses).folder,temp_an(ses).name,'beh.mat'));    
%         nCell = numel(cells.dFF);
%                
%         %% Calculate trial-averaged dF/F and choice selectivity       
%         psth_left = cell(1,nCell);
%         psth_right = cell(1,nCell);
%         psth_left_err = cell(1,nCell);
%         psth_right_err = cell(1,nCell);
%         
%         choicesel_hit = cell(1,nCell);
%         choicesel_err = cell(1,nCell);
%         
%         for j=1:nCell
%             if ~isempty(cells.dFF{j})
%             % calculate PSTH, i.e. trial-averaged dF/F
%             fieldname={'sound','upsweep','left','hit'}; trialMask = getMask(trials,fieldname);
%             psth_left{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
%             % signal is average dff for all events
%             % nEvents = number of all selected trials/events
%             fieldname={'sound','downsweep','right','hit'}; trialMask = getMask(trials,fieldname);
%             psth_right{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
%             
%             fieldname={'sound','downsweep','left','err'}; trialMask = getMask(trials,fieldname);
%             psth_left_err{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
%             fieldname={'sound','upsweep','right','err'}; trialMask = getMask(trials,fieldname);
%             psth_right_err{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
%                   
%             % choicesel_XLeft{j} = calc_selectivity(psth_left{j},psth_left_err{j});        
%             % calculate significancy of choice selectivity
%             %??
% %             choicesel_hit_sig{j} = 0/1;
% %             choicesel_err_sig{j} = 0/1;
%             end
%         end
%         
%         %% plot average choice selectivity for each session
%         subplot(2,1,1); hold on
%         val_sig = nan(nCell, numel(choicesel_hit{1}.signal));
%         for j=1:nCell
%             if ~isempty(choicesel_hit{j})
%             val_sig(j,:) =choicesel_hit{j}.signal;
%             end
%         end
%         t = choicesel_hit{1}.t;
%         plot(t,nanmean(val_sig), 'color', ccmap(ses)); 
%         title ('Hit Trials')
%         box off
%         ylabel ('Time (seconds)')
%         xlabel ('Choice selectivity')
%         
%         subplot(2,1,2); hold on
%         val_sig = nan(nCell, numel(choicesel_hit{1}.signal));
%         for j=1:nCell
%             if ~isempty(choicesel_hit{j})
%             val_sig(j,:) =choicesel_err{j}.signal;
%             end
%         end
%         t = choicesel_hit{1}.t;
%         plot(t,nanmean(val_sig), 'color', ccmap(ses)); 
%         title ('Error Trials')
%         box off
%         ylabel ('Time (seconds)')
%         xlabel ('Choice selectivity')
%         
%     end  
% end
% 
% %% Plot cue-aligned dF/F for a representative cell
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
