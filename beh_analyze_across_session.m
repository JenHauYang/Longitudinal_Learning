% % behavior_longitudinal %%
% OBJECTIVE: analyze behavioral changes across session

% figures save path
fpath='D:\JenHau\siniscalchi2019\Learning\longitudinal\figures';
analysisFilePath = 'D:\JenHau\siniscalchi2019\Learning\analysis';
setup_figprop;

%%  Correct rate analysis
animalList = [{'M52'};{'M53'};{'M54'};{'M55'};{'M56'}];
% ccmap =['k','b','r','g','y']; 
% plot(beh,temp_fig, 's-','color',ccmap(animalID))

%% Analysis for performance
for animalID = 1: numel(animalList)
    curr_animal = animalList{animalID};
    temp_an = dir(fullfile(analysisFilePath,['*',curr_animal,'*']));

    correct_rate = nan(1,size(temp_an,1));
    correct_left = nan(1,size(temp_an,1));
    correct_right = nan(1,size(temp_an,1));
    rt_all = nan(1,size(temp_an,1));
    rt_left = nan(1,size(temp_an,1));
    rt_right = nan(1,size(temp_an,1));
    pre_lick = nan(1,size(temp_an,1));
    lick_bias = nan(1,size(temp_an,1));
    dprime = nan(1,size(temp_an,1));
    wsls_rate = nan(1,size(temp_an,1));
    ws_rate = nan(1,size(temp_an,1));
    ls_rate = nan(1,size(temp_an,1));
    bias = nan(1,size(temp_an,1));
 
    % Get the mean values for each session    
    for ses=1: size(temp_an,1)
        temp = load(fullfile(temp_an(ses).folder,temp_an(ses).name,'beh.mat'));
        
        % calculate correct rate
        output.nTrialsPerformed=sum(~temp.trials.miss);
        corrL = (temp.trials.left & temp.trials.hit);
        corrR = (temp.trials.right & temp.trials.hit);
        correct_rate(ses) = (sum(corrL) + sum(corrR))/output.nTrialsPerformed;
        correct_left(ses) = (sum(corrL)) / sum(temp.trials.upsweep & ~temp.trials.miss); %%sum(temp.trials.upsweep_1 & ~temp.trials.miss); %% ??
        correct_right(ses) = (sum(corrR)) / sum(temp.trials.downsweep & ~temp.trials.miss); 
        
        % calculate RT-  hit & RT vs err& RT // L& hit & RT vs R &hit & RT  
        temp_RT = nan(1, numel(temp.trialData.cueTimes));
        temp_rtL = nan(1, numel(temp.trialData.cueTimes));
        temp_rtR = nan(1, numel(temp.trialData.cueTimes));
        rtL = (temp.trials.left_reaction & temp.trials.hit);
        rtR = (temp.trials.right_reaction & temp.trials.hit);
        
        for k = 1:numel(temp.trialData.cueTimes)
            if ~isnan(temp.trialData.responseTimes(k))
                temp_RT(k) = temp.trialData.responseTimes(k) - temp.trialData.cueTimes(k); % cue presented but mouse did not response
                temp_rtL(k) = rtL(k) * temp_RT(k);
                temp_rtR(k) = rtR(k) * temp_RT(k);
            end
        end
        
        rt_all(ses) = nanmean(temp_RT);
        temp_rtL(temp_rtL == 0) = NaN;
        rt_left(ses) = nanmean(temp_rtL);
        temp_rtR(temp_rtR == 0) = NaN;
        rt_right(ses) = nanmean(temp_rtR);
        
        % calculate link rate
        pre_lick(ses) = nansum(temp.trialData.preCueLickRate);
        lick_bias(ses) = sum(temp.trialData.numLeftLick) / (sum(temp.trialData.numLeftLick) + sum(temp.trialData.numRightLick));
        
        % calculate dprime
        temp_left = correct_left(ses);
        temp_right = correct_right(ses);
        
        if temp_left == 1
            temp_left = 0.99;
        elseif temp_left == 0
            temp_left = 0.01;
        end
        if temp_right == 1
            temp_right = 0.99;
        elseif temp_right == 0
            temp_right = 0.01;
        end
        dprime(ses) = norminv(temp_left) - norminv(1-temp_right);
        
        % calculate WSLS, WS, LS
        ws = temp.trials.hit_1 & ((temp.trials.left_1 & temp.trials.left) | (temp.trials.right_1 & temp.trials.right)); %if they stay after winning prior trial
        ws_denom = temp.trials.hit_1 & (temp.trials.left | temp.trials.right); %only count if they choose on the current trial
        ls = temp.trials.err_1 & ((temp.trials.left_1 & temp.trials.right) | (temp.trials.right_1 & temp.trials.left)); %if they switch after losing prior trial
        ls_denom = temp.trials.err_1 & (temp.trials.left | temp.trials.right); %only count if they choose on the current trial
        
        ws_rate(ses) = sum(ws)./sum(ws_denom); % fraction win-stay out of all wins
        ls_rate(ses) = sum(ls)./sum(ls_denom); % fraction win-stay out of all wins
        wsls_rate(ses) = (sum(ws)+sum(ls))./(sum(ws_denom)+sum(ls_denom)); %fraction win-stay-lose-switch out of all trials
        
        % calculate bias
        bias(ses) = sum(temp.trials.left) / (sum(temp.trials.left) + sum(temp.trials.right));
    end
    
    % Get output for each animal to take average across animals 
    Mouse(animalID).number = curr_animal;
    Mouse(animalID).correct_rate = correct_rate;
    Mouse(animalID).correct_left = correct_left;
    Mouse(animalID).correct_right = correct_right;
    Mouse(animalID).rt_all = rt_all;
    Mouse(animalID).rt_left = rt_left;
    Mouse(animalID).rt_right = rt_right;
    Mouse(animalID).pre_lick = normalize(pre_lick);
    Mouse(animalID).lick_bias = lick_bias;
    Mouse(animalID).dprime = dprime;
    Mouse(animalID).wsls_rate = wsls_rate;
    Mouse(animalID).ws_rate = ws_rate;
    Mouse(animalID).ls_rate = ls_rate;
    Mouse(animalID).bias = bias;
    
end  

%% Get the data for each session & for each mice
% in the same loop calculate 13 variables for all mices BIG MOUSE structure
mice_beh = struct2cell(Mouse);
correct_rate = cell2mat(squeeze(mice_beh(2,:,:)));
correct_left = cell2mat(squeeze(mice_beh(3,:,:)));
correct_right = cell2mat(squeeze(mice_beh(4,:,:)));
rt_all = cell2mat(squeeze(mice_beh(5,:,:)));
rt_left = cell2mat(squeeze(mice_beh(6,:,:)));
rt_right = cell2mat(squeeze(mice_beh(7,:,:)));
pre_lick = cell2mat(squeeze(mice_beh(8,:,:)));
lick_bias = cell2mat(squeeze(mice_beh(9,:,:)));
dprime = cell2mat(squeeze(mice_beh(10,:,:)));
wsls_rate = cell2mat(squeeze(mice_beh(11,:,:)));
ws_rate = cell2mat(squeeze(mice_beh(12,:,:)));
ls_rate = cell2mat(squeeze(mice_beh(13,:,:)));
bias = cell2mat(squeeze(mice_beh(14,:,:)));

beh = [1:size(temp_an,1)]; % for plotting, x axis 

%% plot for performance: corrrect rate, correct left, correct right,dprime, bias, false alarm
h = figure;
for k=1:5
    subplot(2,3,k);hold on
    if k==1
        temp_fig = correct_rate';
        mean_temp_fig = nanmean(temp_fig,2); 
        tName= 'Correct rate';
        yrange = [0 1];
    elseif k==2 
        temp_fig = correct_left';
        mean_temp_fig = nanmean(temp_fig,2); 
        tName ='Correct left';
        yrange = [0 1];
    elseif k==3
        temp_fig = correct_right';
        mean_temp_fig = nanmean(temp_fig,2); 
        tName ='Correct right';
        yrange = [0 1];
    elseif k==4
        temp_fig = dprime';
        mean_temp_fig = nanmean(temp_fig,2); 
        tName ='d prime';
        yrange = [-1 3];
    elseif k==5
        temp_fig = bias';
        mean_temp_fig = nanmean(temp_fig,2); 
        tName ='Bias';
        yrange = [0.25 0.75];
    end
    plot(beh,temp_fig, '.:','color',[0.5 0.5 0.5]); hold on
    plot(beh, mean_temp_fig,'s-','color','k')
    xlim([0 size(temp_an,1)+1]); box off
    xticks([1:5])
    ylim(yrange)
    xlabel('Session')
    ylabel(tName)   

end
saveas(figure(h),fullfile(fpath,['beh_all_animal_performance.fig']),'fig')
saveas(figure(h),fullfile(fpath,['beh_all_animal_choice-beh.png']),'png')

%% plot outcome: RT, RT left, RT right, pre-cue licks, lick bias
h = figure;
for k=1:5
    subplot(2,3,k);hold on
    if k==1
        temp_fig = rt_all';
        mean_temp_fig = nanmean(temp_fig,2); 
        tName ='RT all';
        yrange = [0 1];
    elseif k==2
        temp_fig = rt_left';
        mean_temp_fig = nanmean(temp_fig,2); 
        tName ='RT left';
        yrange = [0 1];
    elseif k==3
        temp_fig = rt_right';
        mean_temp_fig = nanmean(temp_fig,2); 
        tName ='RT right'; 
        yrange = [0 1];
    elseif k==4
        temp_fig = pre_lick';
        mean_temp_fig = nanmean(temp_fig,2); 
        tName ='Pre-cue licks'; 
        yrange = [-2 2];    
    elseif k==5
        temp_fig = lick_bias';
        mean_temp_fig = nanmean(temp_fig,2); 
        tName ='Lick bias'; 
        yrange = [0.25 0.75];
    end
    plot(beh,temp_fig, '.:','color',[0.5 0.5 0.5]); hold on
    plot(beh, mean_temp_fig,'s-','color','k')
    xlim([0 size(temp_an,1)+1]); box off
    xticks([1:5])
    ylim(yrange)
    xlabel('Session')
    ylabel(tName)   
    
end
saveas(figure(h),fullfile(fpath,['beh_all_animal_outcome.fig']),'fig')
saveas(figure(h),fullfile(fpath,['beh_all_animal_choice-beh.png']),'png')

%% plot of choice behaviour: wsls, ws, ls
h = figure;
for k=1:3
    subplot(2,3,k);hold on
    if k==1
        temp_fig = wsls_rate';
        mean_temp_fig = nanmean(temp_fig,2); 
        tName ='WSLS rate';
        yrange = [0.3 0.7];
    elseif k==2
        temp_fig = ws_rate';
        mean_temp_fig = nanmean(temp_fig,2);
        tName ='Win-stay rate';
        yrange = [0.3 0.7];
    elseif k==3
        temp_fig = ls_rate';
        mean_temp_fig = nanmean(temp_fig,2);
        tName ='Lose-shift rate';
        yrange = [0.3 0.7];
    end
    plot(beh,temp_fig, '.:','color',[0.5 0.5 0.5]); hold on
    plot(beh, mean_temp_fig,'s-','color','k')
    xlim([0 size(temp_an,1)+1]); box off
    xticks([1:5])
    ylim(yrange)
    xlabel('Session')
    ylabel(tName)    
end
saveas(figure(h),fullfile(fpath,['beh_all_animal_choice-beh.fig']),'fig')
saveas(figure(h),fullfile(fpath,['beh_all_animal_choice-beh.png']),'png')

%% Annalysis of variance - ls-rate
% data = dprime(:);
% Y = [ones(5,5).*[1:5]];  % session,
% Y = Y(:)
% [p, atab] = anova1(data, Y, 'varnames',{'Session'});
