 function createPSTHsforEachCellforEachAnimal
% this function creates PSTHs for every cell for each animals, takes time
% to finish (~2hrs)
 
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
        
        fname_save = fullfile(temp_an(ses).folder,temp_an(ses).name, 'choice selectivity');
        if ~exist(fname_save,'dir')
            mkdir(fname_save);
        end
        
        params=[];
        params.trigTime = trialData.cueTimes;
        params.xtitle = 'Time from stimulus (s)';
        params.window = [-2:0.25:6.5];
        params.minNumTrial = 5;
        
        %% PSTH for trials in the session for all cells - create PSTH
        psth_panel =[];
        for j = 1: nCell
            if ~isempty(cells.dFF{j})
                for k=1:2 % for panels
                    fieldname=[];
                    if k==1 %panel 1 - Hits
                        fieldname{1}={'sound','upsweep','left','hit'}; col{1}='r'; linstyle{1}='-';
                        fieldname{2}={'sound','downsweep','right','hit'}; col{2}='b'; linstyle{2}='-';
                    elseif k==2 %panel 2 - Errs
                        fieldname{1}={'sound','downsweep','left','err'}; col{1}='r'; linstyle{1}=':';
                        fieldname{2}={'sound','upsweep','right','err'}; col{2}='b'; linstyle{2}=':';
                    end
                    for kk=1:numel(fieldname)
                        trialMask = getMask(trials,fieldname{kk});
                        psth_panel(k).sig{kk} = get_psth( cells.dFF{j}, cells.t, params.trigTime(trialMask), strjoin(fieldname{kk}), params );
                        psth_panel(k).col{kk} = col{kk};
                        psth_panel(k).linstyle{kk} = linstyle{kk};
                    end
                end
                tlabel = ['Cell ' int2str(j)];
                plot_psth(psth_panel,tlabel,params.xtitle);
                print(gcf,'-dpng',fname_save);
                fname_save = fullfile(temp_an(ses).folder,temp_an(ses).name, 'choice selectivity' ,['cell' int2str(j) '-choice']);
                saveas(gcf, fname_save, 'fig');
                close all;
            end
        end   
        
    end
end
toc

% To visualize the mean PSTHs 
% mice_psthAll = struct2cell(dat);
% correct_rate = cell2mat(squeeze(mice_beh(2,:,:)));
% correct_left = cell2mat(squeeze(mice_beh(3,:,:)));
