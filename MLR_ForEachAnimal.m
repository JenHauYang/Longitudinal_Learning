function MLR_ForEachAnimal

%% MLR: 5 sessions for each of 5 animals
% 25 figures: 5 sessions x 5 animals
% calculate and save ~ 1 hour
% cell number in animal 1 = 220; 2 = 266; 3 = 235; 4 = 184; 5 = 199

%% to calculate and save
root_path = 'D:\JenHau\siniscalchi2019\Learning\longitudinal';
data_path = fullfile(root_path, 'data','output');
fig_path = fullfile(root_path, 'figures');
dffFilePath = 'D:\JenHau\siniscalchi2019\Learning\analysis';
animalList = [{'M52'};{'M53'};{'M54'};{'M55'};{'M56'}];

tic
for animalID = 1: numel(animalList)
    curr_animal = animalList{animalID};
    temp_an = dir(fullfile(dffFilePath,['*',curr_animal,'*']));
    
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
        trialMask = getAnyMask(trials,fieldname);
        
        for j=1:numel(cells.dFF)
            reg_cr{j}=linear_regr( cells.dFF{j}, cells.t, [params.choiceEvent params.outcomeEvent], params.trigTime, trialMask, params );
        end
        
        MLRforEach(animalID).session{ses} = reg_cr; 
    end
end
save ( fullfile(data_path,'data_MLRforEach'),'MLRforEach')
toc

%% to plot
setup_figprop;
root_path = 'D:\JenHau\siniscalchi2019\Learning\longitudinal';
fig_path = fullfile(root_path, 'figures');
animalList = [{'M52'};{'M53'};{'M54'};{'M55'};{'M56'}];

load('D:\JenHau\siniscalchi2019\Learning\longitudinal\data\output\data_MLRforEach.mat')

params.xtitle = 'Time from stimulus (s)';
params.pvalThresh = 0.01; 
tlabel={'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'};

for animalID = 1: numel(animalList)
    curr_animal = animalList{animalID};
    for ses=1: length(MLRforEach(1).session)        
        reg_cr = MLRforEach(animalID).session{ses};
        plot_regr(reg_cr,params.pvalThresh,tlabel,params.xtitle);
        print(gcf,'-dpng',fig_path);
        saveas(gcf, fullfile(fig_path,[(curr_animal),'session',num2str(ses)]),'fig');
        saveas(gcf, fullfile(fig_path,[(curr_animal),'session',num2str(ses)]),'png');
        close all;
    end
end

end