clear all;
close all;

% SET ---------------------------------------------------------------------

% Add directories

addpath(genpath('/Users/admin/Documents/MATLAB/'));
cd('/Users/admin/Documents/MATLAB/BackwardsLetters/Study_2/Exp1_Random_0and180degs/Data');

% Task parameters

letterArray = [char(65:66) char(68:85) char(87:90)];      % A to Z
nLeadTrailItems = 6;            % Nuumber of items in stream at beginning and end with no target
nConditions = 2;
nStreams = 2;
nTrials = 100;                   % Per participant, per condition, per session
nSessions = 1;
nBlocks = 4;
itemRate = 8;
% !!!!!!!!!!!!!!

% Participant details

%allParticipants =  {'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AO', 'AP', 'AQ', 'AR', 'AU', 'AW', 'AX', 'BC', 'BD'};
allParticipants = {'AA', 'AB', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AJ', 'AK', 'AL', 'AM', 'AN', 'AO', 'AP', 'AQ', 'AR', 'AS', 'AT', 'AU', 'AV', 'AW', 'AX', 'AY', 'AZ', 'BA', 'BB', 'BC', 'BD', 'BE', 'BF', 'BG', 'BH', 'BI', 'BJ', 'BK', 'BL', 'BM' ,'BN' ,'BO' ,'BP' ,'BQ' ,'BR', 'BS'};

% Model fitting parameters

nFreeParameters = 3;
pdf_normmixture = @pdf_Mixture_Exp1;
pdf_uniformonly = @pdf_Uniform_Exp1;
nReplicates = 100;
pCrit = .05;
smallNonZeroNumber = 10^-10;
fitMaxIter = 10^4;
fitMaxFunEvals = 10^4;


% CALCULATE ---------------------------------------------------------------

% Task parameters
totalTrials = (nTrials*nConditions)/2;
nLetters = length(letterArray); % Number of possible letters
maxError = nLetters-(nLeadTrailItems+1);
errorValues = -maxError:1:maxError;
offsetFactor = maxError+1;
nErrorValues = numel(errorValues);

% Participant details

nParticipants = numel(allParticipants);


% COMPILE -----------------------------------------------------------------

% Set up data structures

compiledErrors = NaN(nConditions,nParticipants,nSessions,nTrials,nStreams);
compiledTargets = NaN(nConditions,nParticipants,nSessions,nTrials,nStreams);
compiledResponses = NaN(nConditions,nParticipants,nSessions,nTrials,nStreams);
compiledResponsePositions = NaN(nConditions,nParticipants,nSessions,nTrials,nStreams);

allRates = NaN(1,nConditions);

% Find relevant files

fileList = what;
fileList = fileList.mat;


for thisParticipant = 1:nParticipants
    
        thisPrefix = allParticipants{thisParticipant};
        loadFiles = find(strncmp(thisPrefix,fileList,numel(thisPrefix)))';

        % Load and compile

        for thisFile = 1:length(loadFiles)

            load(fileList{loadFiles(thisFile)});
            
            trialCounter = zeros(1,nConditions);

            
            for thisTrial = 1:totalTrials

                    % Calculate

                    thisCondition = allConditions(thisTrial);

                    trialCounter(thisCondition) = trialCounter(thisCondition) + 1;

                    targetPosition = allTargets(thisTrial,:);
                    responseLetters = allResponses(thisTrial,:);

                    targetLetters = NaN(1,2);
                    responsePosition = NaN(1,2);

                                for thisStream = 1:nStreams

                                    if ~isnan(targetPosition(thisStream))

                                        targetLetters(thisStream) = allLetterOrder(thisTrial,thisStream,targetPosition(thisStream));

                                    end

                                    if ~isnan(responseLetters(thisStream))

                                        responsePosition(thisStream) = find(squeeze(allLetterOrder(thisTrial,thisStream,:))==responseLetters(thisStream));

                                    end

                                end

                                positionError = responsePosition-targetPosition;
                % Store

                 compiledErrors(thisCondition,thisParticipant,thisFile,trialCounter(thisCondition),:) = positionError;
                 compiledTargets(thisCondition,thisParticipant,thisFile,trialCounter(thisCondition),:) = targetLetters;
                 compiledResponses(thisCondition,thisParticipant,thisFile,trialCounter(thisCondition),:) = responseLetters;
                 compiledResponsePositions(thisCondition,thisParticipant,thisFile,trialCounter(thisCondition),:) = responsePosition;
                 
            end
            
            compiledErrors; %can check errors

        end

%        allRates(thisCondition) = itemRate;

end
fprintf('\n\nCompile complete.\n\n');

% MODEL -------------------------------------------------------------------

% Separated by participant

% Build data structures
allEstimates_byParticipant = NaN(nConditions,nStreams,nParticipants,nFreeParameters);
allLowerBounds_byParticipant = NaN(nConditions,nStreams,nParticipants,nFreeParameters);
allUpperBounds_byParticipant = NaN(nConditions,nStreams,nParticipants,nFreeParameters);

% Set options
options = statset('MaxIter', fitMaxIter, 'MaxFunEvals', fitMaxFunEvals, 'Display', 'off');

for thisCondition = 1:nConditions
    
    fprintf('\n\nCondition %d of %d...', thisCondition, nConditions);
    
    for thisStream = 1:nStreams
        
        fprintf('\nStream %d of %d...', thisStream, nStreams);
        
        for thisParticipant = 1:nParticipants
            
            fprintf('\nParticipant %d of % d', thisParticipant, nParticipants);
       
            minNegLogLikelihood = inf;

             if thisCondition==1
                % Dual-stream, all trials
                theseErrors = squeeze(compiledErrors(1,thisParticipant,:,:,thisStream));
                theseErrors = theseErrors(:);
                theseErrors = theseErrors(~isnan(theseErrors));
                              %  figure;hist(theseErrors,-18:18);
             elseif thisCondition==2
                % Dual-stream, all trials
                theseErrors = squeeze(compiledErrors(2,thisParticipant,:,:,thisStream));
                theseErrors = theseErrors(:);
                theseErrors = theseErrors(~isnan(theseErrors));
                            %    figure;hist(theseErrors,-18:18);
                
             end
            
             % Compute negative log likelihood for uniform distribution
        
            uniformNegLogLikelihood = -sum(log(pdf_uniformonly(theseErrors,1)));

            warning('off', 'stats:mlecov:NonPosDefHessian');

            fprintf('      ');
            
            for thisReplicate = 1:nReplicates
                
                fprintf('\b\b\b\b\b\bR%02d...', thisReplicate);

                pGuess = rand;
                muGuess = ((3*rand)-1.5);
                sigmaGuess = 2*rand;
                parameterGuess = [pGuess muGuess sigmaGuess];
                parameterLowerBound = [0 -4 smallNonZeroNumber];
                parameterUpperBound = [1 4 5];
                
                [currentEstimates, currentCIs] = mle(theseErrors, 'pdf', pdf_normmixture, 'start', parameterGuess, 'lower', parameterLowerBound, 'upper', parameterUpperBound, 'options', options);

                % Compute negative log likelihood
                thisNegLogLikelihood = -sum(log(pdf_normmixture(theseErrors,currentEstimates(1),currentEstimates(2),currentEstimates(3))));

                if minNegLogLikelihood > thisNegLogLikelihood
                    minNegLogLikelihood = thisNegLogLikelihood;
                    bestEstimates = currentEstimates;
                    bestEstimateCIs = currentCIs;
                end

            end

            warning('on', 'stats:mlecov:NonPosDefHessian');

            
            % Test for a significant difference in log likelihoods
            [h,pValue,stat,cValue] = lratiotest(-minNegLogLikelihood,-uniformNegLogLikelihood,nFreeParameters,pCrit);
        
            if h==0
            
                % Null model not rejected; use uniform only

                allEstimates_byParticipant(thisCondition,thisStream,thisParticipant,:) = [0 NaN NaN];
                allLowerBounds_byParticipant(thisCondition,thisStream,thisParticipant,:) = [0 NaN NaN];
                allUpperBounds_byParticipant(thisCondition,thisStream,thisParticipant,:) = [0 NaN NaN];

            else

                % Use mixture

                allEstimates_byParticipant(thisCondition,thisStream,thisParticipant,:) = bestEstimates;
                allLowerBounds_byParticipant(thisCondition,thisStream,thisParticipant,:) = bestEstimateCIs(1,:);
                allUpperBounds_byParticipant(thisCondition,thisStream,thisParticipant,:) = bestEstimateCIs(2,:);
                sca
            end

        end
        
    end
    
end


cd('/Users/admin/Documents/MATLAB/BackwardsLetters/Study_2/Exp1_Random_0and180degs/Data');
save('Exp1_Random_0and180degs_ModelRotated.mat','allEstimates_byParticipant','allLowerBounds_byParticipant','allUpperBounds_byParticipant', 'compiledErrors', 'positionError');

fprintf('\n\n 0 degrees (efficacy)');
squeeze(allEstimates_byParticipant(1,:,:,1))'
fprintf('180 degrees (efficacy)');
squeeze(allEstimates_byParticipant(2,:,:,1))'


fprintf('0 degrees (latency)');
squeeze(allEstimates_byParticipant(1,:,:,2))'
fprintf('180 degrees (latency)');
squeeze(allEstimates_byParticipant(2,:,:,2))'

fprintf('0 degrees (precision)');
squeeze(allEstimates_byParticipant(1,:,:,3))'
fprintf('180 degrees (precision)');
squeeze(allEstimates_byParticipant(2,:,:,3))'
