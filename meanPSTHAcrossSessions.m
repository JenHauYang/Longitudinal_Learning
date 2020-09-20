function meanPSTHAcrossSessions
% this functions calculate mean PSTHs across sessions per animal, takes about ~4
% minutes to finish. Plot means returns snake plot (heatmap).

dffFilePath = 'D:\JenHau\siniscalchi2019\Learning\analysis';
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
    fields = [ {'psth_left'},{'psth_right'},{'psth_left_err'},{'psth_right_err'}];

    for ses=1: size(temp_an,1)
        load(fullfile(temp_an(ses).folder,temp_an(ses).name,'dff.mat'));
        load(fullfile(temp_an(ses).folder,temp_an(ses).name,'beh.mat'));
        nCell = numel(cells.dFF);
        params.trigTime = trialData.cueTimes;
        
        %% Calculate PSTH for all cells
        psth_left = [];
        psth_right = [];
        psth_left_err = [];
        psth_right_err = [];
        
        cIndx =1;
        for j=1:nCell
            if ~isempty(cells.dFF{j})
                % calculate PSTH, i.e. trial-averaged dF/F
                fieldname={'sound','upsweep','left','hit'}; trialMask = getMask(trials,fieldname);
                psth_left{cIndx} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
                % signal is average dff for all events
                % nEvents = number of all selected trials/events
                fieldname={'sound','downsweep','right','hit'}; trialMask = getMask(trials,fieldname);
                psth_right{cIndx} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
                
                fieldname={'sound','downsweep','left','err'}; trialMask = getMask(trials,fieldname);
                psth_left_err{cIndx} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
                fieldname={'sound','upsweep','right','err'}; trialMask = getMask(trials,fieldname);
                psth_right_err{cIndx} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname), params);
                cIndx = cIndx +1; 
            end
        end
        % keep the psth values for each session
        dat(animalID).([fields{1},num2str(ses)]) = psth_left;
        dat(animalID).([fields{2},num2str(ses)]) = psth_right;
        dat(animalID).([fields{3},num2str(ses)]) = psth_left_err;
        dat(animalID).([fields{4},num2str(ses)]) = psth_right_err;
        
    end
end
toc

%% Plot means for each session
fpath='D:\JenHau\siniscalchi2019\Learning\longitudinal\figures\Trial-averaged activity';

% left_hit vs. right_hit
for animalID = 1: numel(animalList)
    curr_animal = animalList{animalID};
    h = figure('Name',(curr_animal),'NumberTitle','off');
    params.xtitle = 'Time from stimulus (s)';
    for ses =1:5
        subplot(2,5,ses)
        params.title = ['Session: ',num2str(ses, '%.2d')];
        temp = dat(animalID).(['psth_left',num2str(ses)]);
        plot_snake(temp,[0 6.5],params);
        
        subplot(2,5,ses+5)
        temp = dat(animalID).(['psth_right',num2str(ses)]);
        plot_snake(temp,[0 6.5],params);
    end
    sgtitle('PSTH left vs. right');
    
    saveas(figure(h),fullfile(fpath,(curr_animal)),'fig')
    saveas(figure(h),fullfile(fpath,(curr_animal)),'png')
    
end

