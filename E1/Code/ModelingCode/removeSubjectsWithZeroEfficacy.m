clear all;
close all;

cd('/Users/admin/Documents/MATLAB/BackwardsLetters/Study_2/Exp1_Random_0and180degs/Data');
load('Exp1_Random_0and180degs_ModelRotated.mat');

allParticipants = {'AA', 'AB', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AJ', 'AK', 'AL', 'AM', 'AN', 'AO', 'AP', 'AQ', 'AR', 'AS', 'AT', 'AU', 'AV', 'AW', 'AX', 'AY', 'AZ', 'BA', 'BB', 'BC', 'BD', 'BE', 'BF', 'BG', 'BH', 'BI', 'BJ', 'BK', 'BL', 'BM' ,'BN' ,'BO' ,'BP' ,'BQ' ,'BR', 'BS'};
participantsToExclude = [5 6 7 13 14 16 17 18 19 20 25 28 34 38 39];

subjectsForAnalysis = allParticipants;
subjectsForAnalysis(participantsToExclude) = [];

E1_allEstMinusExcluded = allEstimates_byParticipant;
E1_allEstMinusExcluded(:,:,participantsToExclude,:) = [];

cd('/Users/admin/Documents/MATLAB/BackwardsLetters/Study_2/Exp1_Random_0and180degs/Data');
save('Exp1_InvertedLtrs_ParamsForAnalysis.mat','E1_allEstMinusExcluded','subjectsForAnalysis');