% master code for longitudunal analysis

%To do for this week:

% 1) Create function for 'data_dffSignalComparison' - run/ save
% 2) Another function - For all cells , create left/right hit signals
% 3) Mean for all cells from one animal + SEM ( 5 figures) 
% 4) Mean for all cells from all animals ( 1 figure)
% 5) Look at peaks magnitude for 5 sesions ( mean +std) - mean for each animal + mean
% of all animal
% 6) Look at the peak locs - ?!!
% 7) Calculate the peak
% 8) run Log reg data structure & save it "stats_mat"
% 9) start writing plotting figure for one animal
% Repository for Jen Hau's code
% Tues: We can spatial localisation changes over time. 

% Check mean choice selectivity for mean of all animals - 5 session
% figures.


%%
% analyze the behaviour 
beh_analyze_across_session
% - add analysis ANOVA or something else? 
% RT good to keep in mind one animal 5th session problem!

%% PSTH  - Average trial activiies
createPSTHsforEachCellforEachAnimal
% creates one figure for each cell, takes ~ 2/3 hours

meanPSTHAcrossSessions
% creates mean PSTHs and plot heatmaps

% latency analysis -( Hit session 1 x 5 sessions)  Hit vs Error 

% 1) change in cell response - dff _real signal
% -  latency 
% -  amplitude
% -  hit vs err representation (choice selectivity) - # of cells, 

% 2) change in population response - # of sig. choice selectivty


                                    
%%  Choice senstivity - across trials across session 
%Check choice selectivity across sessions  - latter session, higher
% selectivity! - copy & paste code snippet from cortex paper
choice_selectivity_across_session
% linear regression to determine choice selective cells, and plot numbers
% of thses cells across sessions
% Change thresholf & look for other contributors (N-1) etc

meanChoiceSelectivityAcrossSessions 
% creates mean choice selectivity and plot heatmaps - (5 sessions

% look at latency analyses in choice selectivity?  plot the same cell response for 5 sesion on
% top of each other
%% Linear regresssion

linearRegressionForChoiceAndOutCome % to calculate the variable

plot_choiceSelectivityDistributin % to plot the figure

