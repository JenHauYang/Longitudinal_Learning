function meanChoiceSelectivitySingleCell
% ~ 6 minutes

% pairwise ttest to find responsiveness for each cell. Focus on left hit, right hit use getMask function), take 1-sec df/f values 
% before and after the onset for all trials, and run ttest2 function do two ttest (1& 2) - so if a cell is significant in 2 ttest 
% - let's call it for a "responsive cell", it means after onset there is a significant change in df/f signal. ( It is responding something)

%% compare df/f signal for [-2 0] 'pre_onset’ to [0 2]  'post_onset’
dffFilePath = 'D:\JenHau\siniscalchi2019\Learning\analysis';
animalList = [{'M52'};{'M53'};{'M54'};{'M55'};{'M56'}];

tic
for animalID = 1: numel(animalList)
    curr_animal = animalList{animalID};
    temp_an = dir(fullfile(dffFilePath,['*',curr_animal,'*']));
    
    for ses=1: size(temp_an,1)
        load(fullfile(temp_an(ses).folder,temp_an(ses).name,'dff.mat'));
        load(fullfile(temp_an(ses).folder,temp_an(ses).name,'beh.mat'));    
        nCell = numel(cells.dFF);
        ttest_left = nan(1,nCell); 
        ttest_right = nan(1,nCell);
                
        params=[];
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from stimulus (s)';
        params.minNumTrial = 5;

        pre_onset_left = [];
        post_onset_left = [];
        pre_onset_right = [];
        post_onset_right = [];
        
        for j = 1: nCell               
            fieldname={'sound','upsweep','left','hit'}; trialMask = getMask(trials,fieldname);
            params.window = [-2:0.25:0];
            pre_onset_left{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
            params.window = [0:0.25:2];
            post_onset_left{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
            
            fieldname={'sound','downsweep','right','hit'}; trialMask = getMask(trials,fieldname);
            params.window = [-2:0.25:0];
            pre_onset_right{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
            params.window = [0:0.25:2];
            post_onset_right{j} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);                             
            
            ttest_left(j) = ttest2(pre_onset_left{j}.signal, post_onset_left{j}.signal);
            ttest_right(j) = ttest2(pre_onset_right{j}.signal, post_onset_right{j}.signal);
        end
        ttest(ses).left = ttest_left; 
        ttest(ses).right = ttest_right;        
    end
    
    % put all array (5 sessions) into a matrix
    t_output_left =  nan(size(temp_an,1),nCell);
    t_output_right =  nan(size(temp_an,1),nCell);
    
    for i = 1:size(temp_an,1)
        t_output_left(i,:) = ttest(i).left;
        t_output_right(i,:) = ttest(i).right;
    end
    
    % determine which cells are selective (sig. ttest in all two conditions) 
    sig_cell_temp = nan(1,nCell);   
    for j = 1:nCell
%         if sum(t_output_left(:,j)) + sum(t_output_right(:,j)) >= 10 %both t-test in all 5 sessions are significant
%             sig_cell_temp(j) = j;
%         end

%         if t_output_left(1,j) + t_output_right(1,j) >= 2 %both t-test in session 1 are significant
%                 sig_cell_temp(j) = j;
%         end        

        if t_output_left(1,j) + t_output_right(1,j) == 0 % in ses 1 - no sig. 
            if t_output_left(2,j) + t_output_right(2,j) >= 2 %sig. in session 2 
                sig_cell_temp(j) = j;
            end
        end
    end
    mouse(animalID).sig_cell = rmmissing(sig_cell_temp);
    mouse(animalID).sig_cell_number = numel(rmmissing(sig_cell_temp));
    
end
toc 

%Get cells significant in left/righ in ses 1

%% plot representative cells
% plot those cells that were significant in all four t-tests in session 1
% ~7 minutes for 15 plots
fpath='D:\JenHau\siniscalchi2019\Learning\longitudinal\figures\Example cells';
save_path = 'D:\JenHau\siniscalchi2019\Learning\longitudinal\data\output';
ccmap =['k','b','r','g','y']; 
setup_figprop;

% load data and calculate variables for plots
clear data
ind =1;
for animalID = 1%:numel(animalList)
    curr_animal = animalList{animalID};
    temp_an = dir(fullfile(dffFilePath,['*',curr_animal,'*']));
    for i = 1:mouse(animalID).sig_cell_number
        j = mouse(animalID).sig_cell(i);
        for ses = 1:size(temp_an,1)
            load(fullfile(temp_an(ses).folder,temp_an(ses).name,'dff.mat'));
            load(fullfile(temp_an(ses).folder,temp_an(ses).name,'beh.mat'));
            
            params.trigTime = trialData.cueTimes;
            % calculate left hit signal
            psth_panel=[];
            fieldname=[];
            fieldname={'sound','upsweep','left','hit'};
            trialMask = getMask(trials,fieldname);
            psth_panel = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params );
            
            val_sig =psth_panel.signal - nanmean(psth_panel.signal (1:7,:)); % -3:0 normalize
            t = psth_panel.t;
            data.leftHit_signal(ind,:) = val_sig';
            data.leftHit_t(ind,:) = t;
            
            % get peak for left-hit signal
            [ pks, locs] =findpeaks(val_sig);
            data.leftHit_peakMagnitude(ind) = max(pks);
            data.leftHit_peakLoc(ind)      = locs(find(pks==max(pks)));
            
            % calculate right hit signal
            psth_panel=[];
            fieldname=[];
            fieldname={'sound','downsweep','right','hit'};
            trialMask = getMask(trials,fieldname);
            psth_panel = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params );
            
            val_sig =psth_panel.signal - nanmean(psth_panel.signal (1:7,:));
            t = psth_panel.t;
            data.rightHit_signal(ind,:) = val_sig';
            data.rightHit_t(ind,:) = t;
            
            % get peak for right-hit signal
            [ pks, locs] =findpeaks(val_sig);
            data.rightHit_peakMagnitude(ind) = max(pks);
            data.rightHit_peakLoc(ind)       = locs(find(pks==max(pks)));
            
            % add other informations
            data.animal(ind) = {curr_animal};
            data.session(ind) = ses;
            data.cellId(ind) = j;
            ind = ind+1;
        end
    end
end
data.note = ('cells which were significant in session 2 but not in session 1');
save ( fullfile(save_path,'data_dffSignalComparison'),'data')

% to plot
params=[];
params.trigTime = trialData.cueTimes;
params.xtitle = 'Time from sound cue (s)';
params.window = [-3:0.5:7];
params.numBootstrapRepeat = 1000;   %number of repeats for bootstrap (for estimating CI)
params.CI = 0.9;  %confidence interval
params.minNumTrial = 0; %only calc PSTH if there is this number of trials
  
ind = 1;
for animalID = 1%:numel(animalList)
    for i = 1:mouse(animalID).sig_cell_number
        j = mouse(animalID).sig_cell(i);
        curr_animal = animalList{animalID};
        
        h = figure('Name',[(curr_animal),' Cell ',num2str(j)],'NumberTitle','off');
        for ses = 1:size(temp_an,1)
            % get wht we need
            cells = animalData(animalID).data(ses).dff.cells;
            trials = animalData(animalID).data(ses).beh.trials;
            trialData =  animalData(animalID).data(ses).beh.trialData; 
            params.trigTime = trialData.cueTimes;
            
            subplot(2,1,1); hold on
            
            plot(t,val_sig,'color',ccmap(ses));
            mean_signal_left(ind,:) = val_sig';
            
            subplot(2,1,2); hold on 
            psth_panel=[];fieldname=[];
            fieldname={'sound','downsweep','right','hit'};
            trialMask = getMask(trials,fieldname);
            psth_panel = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params );                         

            val_sig =psth_panel.signal - nanmean(psth_panel.signal (1:7,:));
            t = psth_panel.t;
            plot(t,val_sig,'color',ccmap(ses));
             mean_signal_right(ind,:) = val_sig';
             
             ind = ind+1;
                
        end
        
         subplot(2,1,1); hold on
         title ('Left hit')
        box off
        xlabel ('Time (seconds)')
        ylabel ('df/f')
        
         subplot(2,1,2); hold on
         title ('Right hit')
        box off
        xlabel ('Time (seconds)')
        ylabel ('df/f')

        legend('Ses 1', 'Ses 2', 'Ses 3', 'Ses 4', 'Ses 5')
        
        saveas(figure(h),fullfile(fpath,[(curr_animal),' Cell ',num2str(j)]),'fig')
        saveas(figure(h),fullfile(fpath,[(curr_animal),' Cell ',num2str(j)]),'png')
        close all
    end    
end
toc

