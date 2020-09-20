function linearRegressionForChoiceAndOutCome
% analyze by MLR, plot fraction of selective cells. ~1 hour

dffFilePath = 'D:\JenHau\siniscalchi2019\Learning\analysis';
fpath='D:\JenHau\siniscalchi2019\Learning\longitudinal\figures\MLR';
animalList = [{'M52'};{'M53'};{'M54'};{'M55'};{'M56'}];

tic
for animalID = 1%: numel(animalList)
    curr_animal = animalList{animalID};
    temp_an = dir(fullfile(dffFilePath,['*',curr_animal,'*']));
%     fields = [ {'choicesel_hit'},{'choicesel_err'},{'choicesel_left'},{'choicesel_right'}];
%     sel_all = nan(1,size(temp_an,1));
    
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
        
        for j=1:5%numel(cells.dFF)
            reg_cr{j}=linear_regr( cells.dFF{j}, cells.t, [params.choiceEvent params.outcomeEvent], params.trigTime, trialMask, params );
        end
        
        temp(animalID).all{ses} = reg_cr; 
%         tlabel={'C(n)','C(n-1)','C(n-2)','R(n)','R(n-1)','R(n-2)','C(n)xR(n)','C(n-1)xR(n-1)','C(n-2)xR(n-2)'};
%         plot_regr(reg_cr,params.pvalThresh,tlabel,params.xtitle);                             
%         print(gcf,'-dpng',fname_save);
%         saveas(gcf, fullfile(fpath,[(curr_animal),'session',num2str(ses)]),'png');
%         close all;
    end
end
toc     
    


% 
% %% if there were task-modulated cells, plot them
% if sum(mod)>0
%     % Make snake plots, choice-selective cells
%     params.xtitle = 'Time from stimulus (s)';
%     plot_snake(psth_left_array(mod),[0 7],params.xtitle);
% 
%     % Choice selectivity plot
%     params.sortParam=[0 6.5];  %sort the cells based on amplitude of selectivity in this time period    
%     params.colorRange=[-0.8 0.8];
% 
%     tlabel = ['Choice selectivity (n=' int2str(sum(mod)) ' choice-selective cells)'];
%     cellOrder = plot_selectivity(choicesel_hit_array(mod),params.sortParam,tlabel,params.xtitle,params.colorRange);
%     print(gcf,'choice-selectivity','-dpng');    %png format
%     saveas(gcf,'choice-selectivity', 'fig');
%     print(gcf,'choice-selectivity','-depsc','-painters');   %eps format
%     
%     %%
%     compTime=[2 4];  %sort the cells based on selectivity in this time period
%     minNumTrial = 1;
%     tlabel = 'Choice selectivity';
%     plot_psth_comp(choicesel_omit_array(mod),choicesel_double_array(mod),compTime,minNumTrial,tlabel);
%     print(gcf,'choice-selectivity_omitdouble_comp','-dpng');    %png format
%     saveas(gcf,'choice-selectivity_omitdouble_comp', 'fig');
%     print(gcf,'choice-selectivity_omitdouble_comp','-depsc','-painters');   %eps format
% 
% end
% 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%     
%     
%     %% ensemble decoding analysis
% runDecode = true;
% if (runDecode)
%         % create array for ensemble activity [time x cell]
%         dFF_ens=[];
%         for j=1:numel(cells.dFF)
%             if ~isempty(cells.dFF{j})
%                 dFF_ens(:,j) = cells.dFF{j};
%             end
%         end
%         
%         %% Decoding CHOICE, linear classifier using every cell in the ensemble
%         
%         % predict choice; dummy-code: left=-1, right=1, miss=NaN
%         params=[];
%         params.trigEvent=NaN(size(trials.left));
%         params.trigEvent(trials.left) = -1;
%         params.trigEvent(trials.right) = 1;
%         % construct linear classifier to decode choice, using dF/F from ensemble, use signals simultaneously recorded
%         params.traintest_fieldname={'hit'};     % use these trials to train classifier and then test (x-fold cross-validation)
% %         params.addtest_fieldname{1}={'doublereward'};    % use the same classifier to additionally test these trials
% %         params.addtest_fieldname{2}={'omitreward'};      % use the same classifier to additionally test these trials
%         params.addtest_fieldname{1}={'err'};    % use the same classifier to additionally test these trials
%         
%         %--- linear classifier
%         params.trigTime = trialData.cueTimes;
%         params.xtitle = 'Time from stimulus (s)';
%         params.window = [-2:0.5:6.5];
%         params.frac = 0.8;      %use this fraction to construct classifier, save the rest for testing, x-fold cross-validation
%         params.numRep = 100;     %number of repeats, should be >=30
%         params.nBack = 1;       %look at decoding of choice for trials n=0 and 1 back
%         
%         lclass_ens_choice=decode_linearclassifier( dFF_ens, cells.t, trials, params );
%         tlabel='Linear classifier';
%         plot_decode(lclass_ens_choice,params.xtitle,tlabel);
%         
% %         %--- same settings, but use random forest classifier
% %         params.numTrees = 100;  %forest with 100 trees
% %         
% %         RF_ens_choice=decode_randomforest( dFF_ens, cells.t, trials, params );
% %         
% %         tlabel='Random forest';
% %         plot_decode(RF_ens_choice,params.xtitle,tlabel);
% %         
% %         %--- same settings, but linear classifier based on leave-one-out sampling
% %         params.frac = [];
% %         params.numRep = [];
% %         lclass_ens_choiceLOO=decode_linearclassifierLOO( dFF_ens, cells.t, trials, params );
% %         tlabel='Linear classifier (LOO)';
% %         plot_decode(lclass_ens_choiceLOO,params.xtitle,tlabel);
% %         
% %         save(fullfile(savematpath,'dff_decode.mat'),...
% %             'lclass_ens_choice','RF_ens_choice','lclass_ens_choiceLOO');
%         
%         %% Decoding CHOICE, adding one cell to ensemble at a time
%         % --- add cell picked randomly
%         numDraw = 100;    %set to low to test the code, for actual analysis should be >=30
%         cellNumList = [1:1:30];  %characterize ensembles with these numbers of cells
%         
%         % predict choice; dummy-code: left=-1, right=1, miss=NaN
%         params=[];
%         params.trigEvent=NaN(size(trials.left));
%         params.trigEvent(trials.left) = -1;
%         params.trigEvent(trials.right) = 1;
%         % construct linear classifier to decode choice, using dF/F from ensemble, use signals simultaneously recorded
%         params.traintest_fieldname={'hit'};     % use these trials to train classifier and then test (x-fold cross-validation)
%         
%         %--- linear/random forest classifier
%         params.trigTime = trialData.cueTimes;
%         params.xtitle = 'Time from stimulus (s)';
%         params.window = [2 4];  %but only look at one time interval
%         params.frac = 0.8;      %use this fraction to construct classifier, save the rest for testing, x-fold cross-validation
%         params.numRep = 1;      %we are making multiple draws, that's where the bootstrap comes from, so will only do 1 5-fold cross-validation for each draw
%         params.numTrees = 100;  %forest with 100 trees
%         
%         lclass_subset_choice = []; RF_subset_choice = [];
%         for j=1:numel(cellNumList)
%             disp(['Ensemble decoding with #cell = ' int2str(cellNumList(j))]);
%             for k=1:numDraw
%                 
%                 idx = randperm(numel(cells.dFF));
%                 dFF_ens_subset = dFF_ens(:,idx(1:cellNumList(j)));  %ensemble constructed from a few randomly drawn cells
%                 
%                 % --- linear classifier
%                 temp = decode_linearclassifier( dFF_ens_subset, cells.t, trials, params );
%                 lclass_subset_choice.corrPred(j,k) = nanmean(temp.corrPred);
%                 lclass_subset_choice.corrPred_randsig(j,k) = nanmean(temp.corrPred_randsig);
%                 lclass_subset_choice.corrPred_scram(j,k) = nanmean(temp.corrPred_scram);
%                 
%                 % --- random forest
%                 temp = decode_randomforest( dFF_ens_subset, cells.t, trials, params );
%                 RF_subset_choice.corrPred(j,k) = nanmean(temp.corrPred);
%                 RF_subset_choice.corrPred_randsig(j,k) = nanmean(temp.corrPred_randsig);
%                 RF_subset_choice.corrPred_scram(j,k) = nanmean(temp.corrPred_scram);
%             end
%             
%             lclass_subset_choice.cellNum(j) = cellNumList(j);
%             RF_subset_choice.cellNum(j) = cellNumList(j);
%         end
%         lclass_subset_choice.fieldname = temp.fieldname;
%         lclass_subset_choice.decode_time = temp.decode_time;
%         lclass_subset_choice.numDraw = numDraw;
%         lclass_subset_choice.numRepeat = temp.numRepeat;
%         RF_subset_choice.fieldname = temp.fieldname;
%         RF_subset_choice.decode_time = temp.decode_time;
%         RF_subset_choice.numDraw = numDraw;
%         RF_subset_choice.numRepeat = temp.numRepeat;
%         
%         tlabel='C(n) - linear classifier';
%         plot_decode_subset(lclass_subset_choice,tlabel);
%         print(gcf,'-dpng',['lclass-subset-choice']);    %png format
%         saveas(gcf, ['lclass-subset-choice'], 'fig');
%         
%         tlabel='C(n) - random forest';
%         plot_decode_subset(RF_subset_choice,tlabel);
%         print(gcf,'-dpng',['RF-subset-choice']);    %png format
%         saveas(gcf, ['RF-subset-choice'], 'fig');
%         
%         save(fullfile(savematpath,'dff_decode.mat'),...
%             'lclass_subset_choice','RF_subset_choice','-append');
%         
%         %% Pairwise correlation, for real and pseudo-ensemble
%         
%         %--- quantify, all cells
%         
%         params=[];        
%         params.numRepeat = 1000;  %number of times to scramble
%         
%         params.trigEvent=NaN(size(trials.left)); %dummy-code: left/single-reward=-1, right/single-reward=1, miss=NaN
%         params.trigEvent(trials.left & trials.hit) = -1;
%         params.trigEvent(trials.right & trials.hit) = 1;
%         
%         params.trigTime = trialData.cueTimes;
%         params.window = [2 4];  %but only look at one time interval
%         
%         pcorr = pairwiseCorr( dFF_ens, cells.t, params );
%         plot_pcorr(pcorr);
%         print(gcf,'pcorr','-dpng');    %png format
%         saveas(gcf,'pcorr', 'fig');
%         
%         %--- quantify, using choice-selective cells
%         
%         pcorr_modCells = pairwiseCorr( dFF_ens(:,mod), cells.t, params );
%         plot_pcorr(pcorr_modCells);
%         print(gcf,'pcorr_modCells','-dpng');    %png format
%         saveas(gcf,'pcorr_modCells', 'fig');
%         
%         save(fullfile(savematpath,'dff_decode.mat'),...
%             'pcorr','pcorr_modCells','-append');
%         
% end
%     
%     %%
% %     clearvars -except i dirs expData expParam runDecode;
% %     close all;

        