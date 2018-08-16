% Dual RSVP stream demo code  
% PTG 11/07/12

Screen('CloseAll')
clear all; 

%Screen('Preference', 'SkipSyncTests', 1)

addpath(genpath('/Users/supalab/Documents/MATLAB/BackwardsLetters/E2/'));

AssertOpenGL;

% Initialise experiment parameters
saveDirectory = '/Users/supalab/Documents/MATLAB/BackwardsLetters/E2/Data/';
userDirectory = '/Users/supalab/Documents/MATLAB/BackwardsLetters/E2/UserData/';
nTrials = 108;            % Per condition upright/inverted - Must be a multiple of 6 as one third are single target, and 1/6th are left/right.
nPracticeTrials = 24;    % Must be a multipe of 12 as one third are single target, 1/6th are left/right, and 1/12th are canonical/inverted.
nSingleTargets = 3;        % Ratio to dual i.e. 1 single every 3 trials
nBlocks = 4;            % 4
nConditions = 2;
nSessions = 1;          % 2
praticeOrder = 1; % 1 = canonical practice first, 2 = inverted practice first.

% Randomise
randomSeed = sum(100*clock);
thisStream = RandStream.create('mt19937ar', 'seed', randomSeed);
RandStream.setGlobalStream(thisStream); 


while true
    participantID = upper(input('Enter two-letter participant initials: ', 's'));
    if length(participantID)==2
        break
    end
end

oldDirectory = cd(userDirectory);
userDataFile = [participantID '_RotatedDualRSVP_Exp2_UserData.mat'];

if ~exist(userDataFile,'file')
   % Make data file
   % Check number of existing records
   allUserDataFiles = what;
   numUserDataFiles = size(allUserDataFiles.mat,1);
   
   % Assign next session
   sessionOrder = randperm(2);
   nextSession = 1;
   totalTime = [];
   
   % Save
   save(userDataFile,'sessionOrder','nextSession','totalTime');
   skipInstructions = 0;
   
else
    
    load(userDataFile);
    skipInstructions = 1;
    
end

cd(oldDirectory);

if nextSession > nSessions
    fprintf('\n\nAll sessions completed!');
    return;
end

% Initialise equipment parameters
runPriority = 1;        % Priority
viewDist = 570;         % Viewing distance in mm
ppm = 3.37; % 1024px/304mm = 3.37  % 2.70  Pixels per mm (GDM, measured 18/7/12)

% Gamma correction
%   y   b   k
% R
% G
% B

% gammaVals = [   3.2173	0.1618  2.6817;
%                 2.4776  0.0096  6.4590;
%                 2.9307  0.0549  2.3531];
%             
% gammaVal_a = 2.005;

gammaVals = 1./[3.2173 2.4776 2.9307];

% Initialise text parameters
letterArray = [char(65:71) char(74:76) char(81:82) char(84:86) char(89)]; %ABCDEFGJKLQRTUVY %[char(65:66) char(68:85) char(87:90)];
canonLetterArray = letterArray;
nLetters = length(letterArray); % Number of possible letters
letterFont = 'Sloan';           % Presentation font - currently MUST be fixed width for drawing accuracy
instructionFont = 'Menlo';
letterTest = 'O';               % Letter for testing letter subtense
letterNull = '-';               % Placeholder letter for response screen
letterSize = 3;                 % Height of letter in degrees of visual angle
letterPoints = 80;              % Startng font size in points (not important, adjusted later)
letterStyle = 0;                % 0 normal, 1 bold
letterSeparation = sqrt(72);    % Distance between letters in degrees of visual angle
letterEccentricity = 6;         % Eccentricity of letters in degrees of visual angle
letterTheta = [pi 0];           % Radial location of letters
letterAlpha = 1.0;              % Opacity of letters

instructionPoints = 16;         % Size of instruction text
textWrap = 50;                  % Number of characters at which to wrap instruction text

%ANGLE
rotateAngle = 180;

% Initialise other spatial parameters
nStreams = 2;                   % Number of letter streams. Currently can't change this parameter.
texturePix = 64;                % Height and width of colour texture, in pixels, is (2*texturePix)+1
textureSD = .6;                 % Standard deviation of Gabor Gaussian in degrees
cueSubtense = 5;                % Subtense of cue in degrees
cueWidth = .1;                  % Width of cue line in degrees
fixationOn = 1;                 % Fixation on (1) or off (0)
fixSubtense = .25;              % Fixation subtense in degrees

% Initialise spatial parameters for response
vSubtenseResponse = 18;         % Vertical subtense of response array in degrees (letters)
hSubtenseResponse = 18;         % Horizontal subtense of response array
responsePoints = 24;      %18   % Font size in points for response letters
okPoints = 14;                  % Font size for 'OK' on response screen
shrinkFactor = .5;              % Factor by which to shrink rect for response collection

% Initialise temporal parameters
simulTargets = 1;               % Simultaneous targets?
nLeadTrailItems = 4;            % Nuumber of items in stream at beginning and end with no target 6
nRepeats = 1;                   % Number of times to repeat unique items in a stream
itemRates = 8;                 % Item rates 10
itemBlankRatio = 2;             % Ratio of item presentation time to inter-item interval
flipProportion = .5;            % Proportion of inter-flip interval prior to flip to request draw
stimResponseWait = .5;          % Number of seconds to wait from offset of stimulus (including blank interval) and onset of response screen
ITI = 2;                        % Number of seconds to wait from response to onset of warning tone
warnWait = 1;                   % Number of seconds to wait after warning tone to first stimulus onset

% Initialise sound parameters
defaultToneLength = .25;        % Length in seconds of default tone
defaultToneFreq = 440;          % Frequency in Hz of default tone
audioSampleRate = 48000;        % Default audio sample rate

% Initialise colour parameters
blackVal = 0;
greyVal = 0;
whiteVal = .75;
fixVal = .75;
cueVal = .75;
letterVal = .75;
dimVal = 0.25;

% Initialise PsychToolbox Pipeline
screenID = max(Screen('Screens'));
PsychImaging('PrepareConfiguration');
%PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'EnablePseudoGrayOutput');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');

[ptbWindow, winRect] = PsychImaging('OpenWindow', screenID, greyVal);
%PsychColorCorrection('SetEncodingGamma', ptbWindow, gammaVals);
[screenWidth, screenHeight] = RectSize(winRect);
screenCentreX = round(screenWidth/2);
screenCentreY = round(screenHeight/2);
flipInterval = Screen('GetFlipInterval', ptbWindow);

% Initialise PsychSound
InitializePsychSound;
ptbAudioPort = PsychPortAudio('Open', [], [], 0, audioSampleRate, 1);

% Instruction text
instructionText1 = ['A rapid stream of letters will appear' ...
    ' on the left and right of the screen.\n\n' ...
    'The two streams will be simultaneous, but the order of ' ...
    'letters will be different.\n\n' ...
    'Within each stream some of the letters will be rotated 180 degrees.\n\n' ...
    'Click the mouse button for more instructions.'];

instructionText2 = ['White circles will appear around one or both streams.\n\n' ...
    'Your task is to report the circled letter/s and their orientation.\n\n' ...
    'The letters may be any of a set of 16 letters from the alphabet.\n\n' ...
    'Click the mouse button for more instructions.'];

instructionText3 = ['Please keep your eyes fixated on the white dot at ' ...
    'the centre of the screen at all times during the letter stream.\n\n' ...
    'Make your response by clicking on the letter/orientation ' ...
    'that corresponds to the circled letter for that stream\n\n' ...
    'If you are unsure, please make your best guess.\n\n' ...
    'Click the mouse button to begin practice.'];

%instructionText4 = 'Click the mouse button to begin practice.';
    
endText = ['The experiment is complete.\n\n' ...
    'Please click the mouse button to end the program.'];


skippedText = 'Please click the mouse button to begin the experiment.';



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Calculate experiment parameters
totalTrials = nTrials*nConditions;
trialsPerBlock = ceil(totalTrials/nBlocks);
stopTrials = 0:trialsPerBlock:totalTrials;
stopTrials(1) = [];
stopTrials(numel(stopTrials)) = [];

% Calculate equipment parameters
mpd = (viewDist/2)*tan(deg2rad(2*letterEccentricity))/letterEccentricity; % Calculate mm per degree of visual angle to the eccentricity of the stimuli
ppd = ppm*mpd;                                              % Calculate pixels per degree of visual angle

% Calculate spatial parameters
cueSubtense_pix = cueSubtense*ppd;                              % Subtense of cue in pixels
cueWidth_pix = cueWidth*ppd;                                    % Width of cue in pixels (pen)
fixSubtense_pix = fixSubtense*ppd;                              % Subtense of fixation marker in pixels
letterHeight_pix = letterSize*ppd;                              % Height of letter in pixels
letterEccentricity_pix = letterEccentricity*ppd;                % Eccentricity of letters in pixels
letterRad = [letterEccentricity_pix letterEccentricity_pix];
[letterX, letterY] = pol2cart(letterTheta, letterRad);
textureSD_pix = textureSD*ppd;                                  % Standard deviation of Gaussian in pixels

% Adjust font size to meet subtense requested
Screen('TextFont', ptbWindow, letterFont);
Screen('TextSize', ptbWindow, letterPoints);
Screen('TextStyle', ptbWindow, letterStyle);

defaultBounds = Screen('TextBounds',ptbWindow,upper(letterTest));
defaultHeight = defaultBounds(4)-defaultBounds(2);
scaleFactor = letterHeight_pix/defaultHeight;
letterPoints = round(letterPoints*scaleFactor);

% Verify
Screen('TextSize', ptbWindow, letterPoints);
Screen('TextFont', ptbWindow, letterFont);
defaultBounds = Screen('TextBounds',ptbWindow,upper(letterTest));
actualLetterHeight = (defaultBounds(4)-defaultBounds(2))/ppd;

% Caluclate sound parameters
defaultTone = sin(linspace(0,2*pi*defaultToneFreq,audioSampleRate*defaultToneLength));

% Calculate colour parameters
halfRange = whiteVal-greyVal;

% Calculate variables for determining order of presentation
nItems = nLetters*nRepeats;                            % Number of items in RSVP stream
nPossTargets = nItems-(2*nLeadTrailItems);

% Calculate spatial parameters for response screen
vSubtenseResponse_pix = vSubtenseResponse*ppd;
hSubtenseResponse_pix = hSubtenseResponse*ppd;

yResponseArrayLetters_Centre = linspace(-vSubtenseResponse_pix/2,vSubtenseResponse_pix/2,nLetters)+screenCentreY;
xResponseArray = linspace(-hSubtenseResponse_pix/2,hSubtenseResponse_pix/2,9)+screenCentreX;

xResponseArrayLetters_Centre = linspace(-hSubtenseResponse_pix/2,hSubtenseResponse_pix/2,nLetters)+screenCentreX;
yResponseArray = linspace(-vSubtenseResponse_pix/2,vSubtenseResponse_pix/2,9)+screenCentreY;

shrinkFactor_pix = texturePix*shrinkFactor;

% Enable alpha blending for typical drawing of masked textures:
%Screen('BlendFunction', ptbWindow, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Create Gaussian aperture
[xGrid,yGrid] = meshgrid(-texturePix:texturePix, -texturePix:texturePix);
gaussianStim = whiteVal*(exp(-(xGrid.^2)/(2*textureSD_pix^2)).*exp(-(yGrid.^2)/(2*textureSD_pix^2)));
%gStim = whiteVal*ones([size(gaussianStim) 4]);
%gStim(:,:,4) = gaussianStim;
%colourTex = Screen('MakeTexture', ptbWindow, gStim);
colourTex = Screen('MakeTexture', ptbWindow, gaussianStim, [], [], 1);

% Put the default tone into the audio buffer
PsychPortAudio('FillBuffer', ptbAudioPort, defaultTone);

% Allocate destination rectatngles
letterRect = Screen('TextBounds', ptbWindow, upper(letterTest));
Screen('TextSize', ptbWindow, responsePoints);
responseLetterRect = Screen('TextBounds', ptbWindow, upper(letterTest));
Screen('TextSize', ptbWindow, letterPoints);
Screen('TextFont', ptbWindow, letterFont);

letterPosX = letterX+screenCentreX;
letterPosY = letterY+screenCentreY;

fixRect = [0 0 fixSubtense_pix fixSubtense_pix];
fixRect = CenterRectOnPoint(fixRect, screenCentreX, screenCentreY);
colourRect = Screen('Rect', colourTex);
cueRect = [0 0 cueSubtense_pix cueSubtense_pix];

cueRects = NaN(2,4);
tempCueX = unique(letterPosX);
tempCueY = unique(letterPosY);

for thisRect = 1:2
   
    cueRects(thisRect,:) = CenterRectOnPoint(cueRect, letterPosX(thisRect), letterPosY(thisRect));

end

cueRectSelect = [1 2];   % Correspondence between cue rects and letter rects

letterPosX = letterPosX-letterRect(3)/2;
letterPosY = letterPosY-letterRect(4)/2;

yResponseArrayLetters = yResponseArrayLetters_Centre-(responseLetterRect(4)/2);
xResponseArrayLetters = xResponseArrayLetters_Centre-(responseLetterRect(3)/2);

validRespRects = zeros(4,4);
validRespRects(1,:) = CenterRectOnPoint([0 0 responseLetterRect(3) vSubtenseResponse_pix+responseLetterRect(4)], xResponseArray(9)+1.5*responsePoints, screenCentreY); %Left-outer (upright)
validRespRects(2,:) = CenterRectOnPoint([0 0 responseLetterRect(3) vSubtenseResponse_pix+responseLetterRect(4)], xResponseArray(9), screenCentreY); %Left-inner (inverted)
validRespRects(3,:) = CenterRectOnPoint([0 0 responseLetterRect(3) vSubtenseResponse_pix+responseLetterRect(4)], xResponseArray(1), screenCentreY); %Right-inner (upright)
validRespRects(4,:) = CenterRectOnPoint([0 0 responseLetterRect(3) vSubtenseResponse_pix+responseLetterRect(4)], xResponseArray(1)-1.5*responsePoints, screenCentreY); %Right-outer (inverted)

xResponseArray = xResponseArray-[responseLetterRect(3) 0 0 letterRect(3) 0 letterRect(3) 0 0 responseLetterRect(3)]./2;
yResponseArray = yResponseArray-[responseLetterRect(4) 0 0 letterRect(4) 0 letterRect(4) 0 0 responseLetterRect(4)]./2;

goRect = CenterRectOnPoint(colourRect, screenCentreX, screenCentreY);

respShrinkPix = [shrinkFactor_pix shrinkFactor_pix -shrinkFactor_pix -shrinkFactor_pix];
goRectResp = goRect+respShrinkPix;

% Condition order

%Trial/condition order
allConditions = [ones(1,nTrials) repmat(2,1,nTrials)]; % Generate list of trial conditions, equal num per condition
numTargets = repmat([2 2 1], 1, totalTrials/nSingleTargets); % One third are single target
targetSide = repmat([0 0 1 0 0 2], 1, nTrials/nSingleTargets);
% Make random
trialRandomiser = randperm(totalTrials);
allConditions = allConditions(trialRandomiser);
numTargets = numTargets(trialRandomiser);
targetSide = targetSide(trialRandomiser);

% Practise trial/condition order
practiceConditions = [ones(1,nPracticeTrials/nConditions) repmat(2,1,nPracticeTrials/nConditions)]; % Generate list of trial conditions, equal num per condition
practiceNumTargets = repmat([2 2 1], 1, nPracticeTrials/nSingleTargets); % One third are single target
practiceTargetSide = repmat([0 0 1 0 0 2], 1, nPracticeTrials/nSingleTargets*2 );
% Make random
practiceRandomiser = randperm(nPracticeTrials);
practiceConditions = practiceConditions(practiceRandomiser);
practiceNumTargets = practiceNumTargets(practiceRandomiser);
practiceTargetSide = practiceTargetSide(practiceRandomiser);

% --------------------------  
    % Calculate temporal parameters
    itemRate = itemRates(1);                 % Item presentation rate (items per second)
    itemBlankDuration = 1.0/itemRate;                                       % Duration in secs of item + blank
    itemDuration = (itemBlankRatio/(itemBlankRatio+1))*itemBlankDuration;   % Duration in secs of item
    blankDuration = itemBlankDuration-itemDuration;    


% Instruction screen
Screen('TextSize', ptbWindow, instructionPoints);
Screen('TextFont', ptbWindow, instructionFont);

%skipInstructions = 1;

if ~skipInstructions

    DrawFormattedText(ptbWindow, instructionText1, 'center', 'center', whiteVal, textWrap);
    Screen('Flip', ptbWindow);

    DrawFormattedText(ptbWindow, instructionText2, 'center', 'center', whiteVal, textWrap);
    WaitSecs(2);
    GetClicks([],0);
    startTotalTime = GetSecs;
    Screen('Flip', ptbWindow);

    DrawFormattedText(ptbWindow, instructionText3, 'center', 'center', whiteVal, textWrap);
    WaitSecs(2);
    GetClicks([],0);
    Screen('Flip', ptbWindow);

%     DrawFormattedText(ptbWindow, instructionText4, 'center', 'center', whiteVal, textWrap);
%     WaitSecs(2);
%     GetClicks([],0);
% 
%     Screen('Flip', ptbWindow);

    WaitSecs(2);
    GetClicks([],0);
    Screen('Flip', ptbWindow);

    Screen('TextSize', ptbWindow, letterPoints);
    Screen('TextFont', ptbWindow, letterFont);
    
    startBlockTime = GetSecs;
    clickTime = startBlockTime;
    HideCursor;
            
    for thisTrial = 1:nPracticeTrials

            Screen('TextSize', ptbWindow, letterPoints);
            Screen('TextFont', ptbWindow, letterFont);
            
            % Set practise block condtion
            
            thisCondition = practiceConditions(thisTrial);
            
            if simulTargets
                thisTarget = repmat(randi(nPossTargets)+nLeadTrailItems,1,2);
            else
                thisTarget = randi(nPossTargets,[1 2])+nLeadTrailItems;
            end

            letterStream = mod([randperm(nItems); randperm(nItems)],nLetters)+1;
            letterOrientaion =  repmat(0:1,1,nLetters/2); % Must be an even number of letters in letterArray
            letterOrientaion = letterOrientaion(randperm(length(letterOrientaion)));
            
            %Check target has correct orientation
            if (thisCondition==1)&&(letterOrientaion(thisTarget(1,1))==1)
                uprightLetters = find(letterOrientaion==0);
                letterOrientaion(uprightLetters(randperm(length(uprightLetters),1))) = 1;
                letterOrientaion(thisTarget(1,1)) = 0;
            elseif (thisCondition==2)&&(letterOrientaion(thisTarget(1,1))==0)                
                invertedLetters = find(letterOrientaion);
                letterOrientaion(invertedLetters(randperm(length(invertedLetters),1))) = 0;
                letterOrientaion(thisTarget(1,1)) = 1;  
            end
            
            letterX = letterPosX;
            letterY = letterPosY;

            theseCueRects = cueRects';

            thisXresponse = repmat([xResponseArray(1) xResponseArray(9)]',1,nLetters);
            thisYresponse = repmat(yResponseArrayLetters,2,1);

            
            % Increase priority

            oldPriority = Priority(runPriority);

            % Sync to VBL at start of animation loop
            Screen('Flip', ptbWindow);

            if fixationOn
                Screen('FillOval', ptbWindow, fixVal, fixRect);
            end

            Screen('FrameOval', ptbWindow, cueVal, theseCueRects(:,1), cueWidth_pix, cueWidth_pix);
            Screen('FrameOval', ptbWindow, cueVal, theseCueRects(:,2), cueWidth_pix, cueWidth_pix);

            % Play the tone
            PsychPortAudio('Start', ptbAudioPort, 1, clickTime + ITI, 0);

            % Display fixation with tone
            warnTime = Screen('Flip', ptbWindow, clickTime + ITI);

            % Remove cue
            if fixationOn
                Screen('FillOval', ptbWindow, fixVal, fixRect);
            end

            Screen('Flip', ptbWindow, warnTime + defaultToneLength - flipProportion*flipInterval);

            % Draw first frame           
            DrawFormattedText(ptbWindow, letterArray(letterStream(1,1)), letterX(1), 'center', letterVal, [], letterOrientaion(1), letterOrientaion(1));
            DrawFormattedText(ptbWindow, letterArray(letterStream(2,1)), letterX(2), 'center', letterVal, [], letterOrientaion(1), letterOrientaion(1));
                
            if fixationOn
                Screen('FillOval', ptbWindow, fixVal, fixRect);
            end

            % Mark drawing ops as finished
            Screen('DrawingFinished', ptbWindow);

            % Done. Draw item.
            startTime = Screen('Flip', ptbWindow, warnTime + warnWait);

            % Blank screen with fixation

            if fixationOn
                Screen('FillOval', ptbWindow, fixVal, fixRect);
            end

            Screen('Flip', ptbWindow, startTime + itemDuration - flipProportion*flipInterval);

            % Start remainder of RSVP loop

            for thisItem = (2:nItems)
                
                % Draw Text
                DrawFormattedText(ptbWindow, letterArray(letterStream(1,thisItem)), letterX(1), 'center', letterVal, [], letterOrientaion(thisItem), letterOrientaion(thisItem));
                DrawFormattedText(ptbWindow, letterArray(letterStream(2,thisItem)), letterX(2), 'center', letterVal, [], letterOrientaion(thisItem), letterOrientaion(thisItem));
                
                if fixationOn
                    Screen('FillOval', ptbWindow, fixVal, fixRect);
                end

                if (thisItem == thisTarget(1)) && (practiceNumTargets(thisTrial) == 2 || practiceTargetSide(thisTrial) == 1)
                    % Draw cue on left
                    Screen('FrameOval', ptbWindow, cueVal, theseCueRects(:,1), cueWidth_pix, cueWidth_pix);
                end

                if thisItem == thisTarget(2) && (practiceNumTargets(thisTrial) == 2 || practiceTargetSide(thisTrial) == 2)
                    % Draw cue on right
                    Screen('FrameOval', ptbWindow, cueVal, theseCueRects(:,2), cueWidth_pix, cueWidth_pix);
                end

                % Mark drawing ops as finished
                Screen('DrawingFinished', ptbWindow);

                % Done. Draw item.
                Screen('Flip', ptbWindow, startTime + (thisItem-1)*itemBlankDuration - flipProportion*flipInterval);

                % Blank screen with fixation

                if fixationOn&&(thisItem<nItems)
                    Screen('FillOval', ptbWindow, fixVal, fixRect);
                end

                Screen('Flip', ptbWindow, startTime + (thisItem-1)*itemBlankDuration + itemDuration - flipProportion*flipInterval);
            end

            % Restore priority
            Priority(oldPriority);

            WaitSecs(stimResponseWait);

            % Draw response screen
            Screen('TextSize', ptbWindow, responsePoints);
            SetMouse(screenCentreX,screenCentreY, screenID);
            ShowCursor('Arrow');

            % Randomly determine response order
            if practiceTargetSide(thisTrial)==0
                responseOrder = randi(2); % 1=L first, 2=R first
            elseif practiceTargetSide(thisTrial)==1
                responseOrder = 1;
            else 
                responseOrder = 2;  
            end
            
            responseTime = NaN(1,2);
            selectedLetters = zeros(nLetters,2);
            selectedOrientation = NaN(1,2);
            rTimeStart = 0;

            for thisResponse = 1:2
                
                if ~(practiceNumTargets(thisTrial)==1 && thisResponse==2) % Don't loop twice if single target
                    
                    if ((thisResponse==1)&&(responseOrder==1))||((thisResponse==2)&&(responseOrder==2))
                        % Left/Top
                        responseSide = 1;
                        colValL = whiteVal;
                        colValR = dimVal;
                        validRects1 = validRespRects(1,:); 
                        validRects2 = validRespRects(2,:);
                    
                    else
                        % Right/Bottom
                        responseSide = 2;
                        colValL = dimVal;
                        colValR = whiteVal;
                        validRects1 = validRespRects(3,:); 
                        validRects2 = validRespRects(4,:);
                    end
                

                goodResponse = 0; %0

                while goodResponse==0

                    % Check that mouse is still within ptbWindow
                    [xMPos,yMPos] = GetMouse(screenID);

                    if xMPos < 10
                        SetMouse(screenCentreX, yMPos, screenID);
                        ShowCursor('Arrow');
                    end

                % Draw attribute matrix
                for thisLetter = 1:nLetters
                    if practiceNumTargets(thisTrial)==2 || practiceTargetSide(thisTrial)==1
                        %Left line-up
                        DrawFormattedText(ptbWindow, letterArray(thisLetter), (thisXresponse(1,thisLetter)-1.5*responsePoints), thisYresponse(1,thisLetter)+responsePoints, colValL);
                        DrawFormattedText(ptbWindow, letterArray(thisLetter), thisXresponse(1,thisLetter), thisYresponse(1,thisLetter)-responsePoints, colValL, [], true, true);
                    end
                    if practiceNumTargets(thisTrial)==2 || practiceTargetSide(thisTrial)==2
                        %Right line-up
                        DrawFormattedText(ptbWindow, letterArray(thisLetter), thisXresponse(2,thisLetter), thisYresponse(2,thisLetter)+responsePoints, colValR);
                        DrawFormattedText(ptbWindow, letterArray(thisLetter), (thisXresponse(2,thisLetter)+1.5*responsePoints), thisYresponse(2,thisLetter)-responsePoints, colValR, [], true, true);
                    end
                end

                    % Draw response letters
                    Screen('TextSize', ptbWindow, letterPoints);
                    Screen('TextFont', ptbWindow, letterFont);

                    if sum(selectedLetters(:,1))
                        %Screen('DrawText', ptbWindow, letterArray(find(selectedLetters(:,1))), letterX(1), letterY(1), colValL, greyVal); %#ok<FNDSB>
                        %flipLetter = letterOrientaion(find(letterStream(1,:)==find(selectedLetters(:,1)))
                        DrawFormattedText(ptbWindow, letterArray(find(selectedLetters(:,1))), letterX(1), 'center', colValL, [], selectedOrientation(1), selectedOrientation(1)); %#ok<FNDSB>
                    elseif practiceTargetSide(thisTrial)~=2
                        DrawFormattedText(ptbWindow, letterNull, letterX(1), 'center', colValL);
                    end

                    if sum(selectedLetters(:,2))
                        %Screen('DrawText', ptbWindow, letterArray(find(selectedLetters(:,2))), letterX(2), letterY(2), colValR, greyVal); %#ok<FNDSB>
                        %flipLetter = letterOrientaion(find(letterStream(2,:)==find(selectedLetters(:,2))));
                        DrawFormattedText(ptbWindow, letterArray(find(selectedLetters(:,2))), letterX(2), 'center', colValR, [], selectedOrientation(2), selectedOrientation(2)); %#ok<FNDSB>
                    elseif practiceTargetSide(thisTrial)~=1
                        DrawFormattedText(ptbWindow, letterNull, letterX(2), 'center', colValR);
                    end

                    Screen('TextSize', ptbWindow, responsePoints);

                    if sum(selectedLetters(:,responseSide)) % If attribute is selected
                        Screen('DrawTexture', ptbWindow, colourTex, [], goRect, [], [], [], [], []);
                        Screen('TextSize', ptbWindow, okPoints);
                        DrawFormattedText(ptbWindow, 'OK', 'center', 'center', blackVal);
                        Screen('TextSize', ptbWindow, responsePoints);
                    else
                        % Possibly draw black GO texture
                        % Screen('DrawTexture', ptbWindow, colourTex, [], goRect, [], [], [], [blackVal blackVal blackVal], []);
                        % DrawFormattedText(ptbWindow, 'OK', 'center', 'center', whiteVal);
                    end

                    if ~ rTimeStart
                        rTimeStart = Screen('Flip', ptbWindow);
                    else
                        Screen('Flip', ptbWindow);
                    end

                    % Get a click
                    [nClicks,xClick,yClick,buttonClick] = GetClicks(ptbWindow,0);
                    clickTime = GetSecs;

                    xClick = screenWidth-xClick;

                    % Check whether the click is in a letter rect
                    if IsInRect(xClick,yClick,validRects1)||IsInRect(xClick,yClick,validRects2)
                        % Check for closest position (Y)
                        [minDist,clickedLetter] = min(abs(yResponseArrayLetters_Centre-yClick));
                        selectedLetters(:,responseSide) = zeros(nLetters,1);
                        selectedLetters(clickedLetter,responseSide) = 1;
                        if IsInRect(xClick,yClick,validRects1)
                            selectedOrientation(responseSide) = 0;
                        else
                            selectedOrientation(responseSide) = 1;
                        end
                    end

                    if IsInRect(xClick,yClick,goRectResp)
                        if sum(selectedLetters(:,responseSide)) % If all attributes are selected
                            responseTime(thisResponse) = clickTime-rTimeStart;
                            goodResponse=1;
                            Screen('Flip', ptbWindow);
                        else
                        % Possible error tone
                        end
                    end
                end
                end
            end


            if thisTrial < nPracticeTrials && thisTrial ~= nPracticeTrials/2
                nextTrialInfoText = ['Click the mouse button to start the next practice trial.\n\n(Trial ' num2str(thisTrial+1) ' of ' num2str(nPracticeTrials) ')'];
                Screen('TextSize', ptbWindow, instructionPoints);
                Screen('TextFont', ptbWindow, instructionFont);
                DrawFormattedText(ptbWindow, nextTrialInfoText, 'center', 'center', whiteVal);
                HideCursor;
                Screen('Flip', ptbWindow);

                % Wait for click to continue
                GetClicks(ptbWindow,0);

            elseif thisTrial == nPracticeTrials/2
                nextTrialInfoText = ['You have completed half the practice trials.'...
                '\n\nIf you have any questions, please ask the experimenter now.\n\n'...        
                'Click the mouse button to start the next practice trial.\n\n(Trial ' num2str(thisTrial+1) ' of ' num2str(nPracticeTrials) ')'];
                Screen('TextSize', ptbWindow, instructionPoints);
                Screen('TextFont', ptbWindow, instructionFont);
                DrawFormattedText(ptbWindow, nextTrialInfoText, 'center', 'center', whiteVal);
                HideCursor;
                Screen('Flip', ptbWindow);

                % Wait for click to continue
                GetClicks(ptbWindow,0);    
                
            else
                nextTrialInfoText = ['Thank you. The practice trials are complete.'...
                '\n\nIf you have any questions, please ask the experimenter now.\n\n'...
                'Please click the mouse button to begin the experiment.'];
                
                Screen('TextSize', ptbWindow, instructionPoints);
                Screen('TextFont', ptbWindow, instructionFont);
                DrawFormattedText(ptbWindow, nextTrialInfoText, 'center', 'center', whiteVal);
            end

    end

else

    DrawFormattedText(ptbWindow, skippedText, 'center', 'center', whiteVal, textWrap);
    startTotalTime = GetSecs;
end

Screen('Flip', ptbWindow);

WaitSecs(2);
GetClicks([],0);
Screen('Flip', ptbWindow);

Screen('TextSize', ptbWindow, letterPoints);
Screen('TextFont', ptbWindow, letterFont);
    
    % --------------------------  
    % Calculate temporal parameters
    itemRate = itemRates(1);                                                % Item presentation rate (items per second)
    itemBlankDuration = 1.0/itemRate;                                       % Duration in secs of item + blank
    itemDuration = (itemBlankRatio/(itemBlankRatio+1))*itemBlankDuration;   % Duration in secs of item
    blankDuration = itemBlankDuration-itemDuration;                         % Duration in secs of blank
    
    % Set up data structures
    allLetterOrder = NaN(totalTrials,2,nItems);
    allLetterOrientations = NaN(totalTrials,nItems);
    
    allTargets = NaN(totalTrials,2);
    allTargetOrientations = NaN(totalTrials,2);
    
    allResponses = NaN(totalTrials,2);
    allResponseOrientations = NaN(totalTrials,2);
    allResponseOrder = NaN(totalTrials,1);
    
    allRTs = NaN(totalTrials,2);
    allPositions = mod(randperm(totalTrials)-1,2)+1;
    
    
    
    % --------------------------
    % BEGIN TRIAL

    startBlockTime = GetSecs;
    clickTime = startBlockTime;
    HideCursor;

    for thisTrial = 1:totalTrials %!!!!!!!!!!!!!
        Screen('TextSize', ptbWindow, letterPoints);
        Screen('TextFont', ptbWindow, letterFont);
        
        letterStream = mod([randperm(nItems); randperm(nItems)],nLetters)+1;
        letterOrientaion =  repmat(0:1,1,nLetters/2); % Must be an even number of letters in letterArray
        letterOrientaion = letterOrientaion(randperm(length(letterOrientaion)));
        thisCondition = allConditions(thisTrial);
        
        if simulTargets
            thisTarget = repmat(randi(nPossTargets)+nLeadTrailItems,1,2);
        else
            thisTarget = randi(nPossTargets,[1 2])+nLeadTrailItems;
        end
        
       %Check target has correct orientation
        if (thisCondition==1)&&(letterOrientaion(thisTarget(1,1))==1)
            uprightLetters = find(letterOrientaion==0);
            letterOrientaion(uprightLetters(randperm(length(uprightLetters),1))) = 1;
            letterOrientaion(thisTarget(1,1)) = 0;
        elseif (thisCondition==2)&&(letterOrientaion(thisTarget(1,1))==0)                
            invertedLetters = find(letterOrientaion);
            letterOrientaion(invertedLetters(randperm(length(invertedLetters),1))) = 0;
            letterOrientaion(thisTarget(1,1)) = 1;  
        end
        
        letterX = letterPosX;
        letterY = letterPosY;

        theseCueRects = cueRects';
        
        thisXresponse = repmat([xResponseArray(1) xResponseArray(9)]',1,nLetters);
        thisYresponse = repmat(yResponseArrayLetters,2,1);

        % Increase priority

        oldPriority = Priority(runPriority);

        % Sync to VBL at start of animation loop
        Screen('Flip', ptbWindow);

        if fixationOn
            Screen('FillOval', ptbWindow, fixVal, fixRect);
        end
        
        Screen('FrameOval', ptbWindow, cueVal, theseCueRects(:,1), cueWidth_pix, cueWidth_pix);
        Screen('FrameOval', ptbWindow, cueVal, theseCueRects(:,2), cueWidth_pix, cueWidth_pix);
        
        % Play the tone
        PsychPortAudio('Start', ptbAudioPort, 1, clickTime + ITI, 0);

        % Display fixation with tone
        warnTime = Screen('Flip', ptbWindow, clickTime + ITI);
        
        % Remove cue
        if fixationOn
            Screen('FillOval', ptbWindow, fixVal, fixRect);
        end

        Screen('Flip', ptbWindow, warnTime + defaultToneLength - flipProportion*flipInterval);

        % Draw first frame
        DrawFormattedText(ptbWindow, letterArray(letterStream(1,1)), letterX(1), 'center', letterVal, [], letterOrientaion(1), letterOrientaion(1));
        DrawFormattedText(ptbWindow, letterArray(letterStream(2,1)), letterX(2), 'center', letterVal, [], letterOrientaion(1), letterOrientaion(1));

        if fixationOn
            Screen('FillOval', ptbWindow, fixVal, fixRect);
        end

        % Mark drawing ops as finished
        Screen('DrawingFinished', ptbWindow);

        % Done. Draw item.
        startTime = Screen('Flip', ptbWindow, warnTime + warnWait);

        % Blank screen with fixation

        if fixationOn
            Screen('FillOval', ptbWindow, fixVal, fixRect);
        end

        Screen('Flip', ptbWindow, startTime + itemDuration - flipProportion*flipInterval);

        % Start remainder of RSVP loop

        for thisItem = (2:nItems)

            % Draw Text
            DrawFormattedText(ptbWindow, letterArray(letterStream(1,thisItem)), letterX(1), 'center', letterVal, [], letterOrientaion(thisItem), letterOrientaion(thisItem));
            DrawFormattedText(ptbWindow, letterArray(letterStream(2,thisItem)), letterX(2), 'center', letterVal, [], letterOrientaion(thisItem), letterOrientaion(thisItem));

            if fixationOn
                Screen('FillOval', ptbWindow, fixVal, fixRect);
            end

            if (thisItem == thisTarget(1)) && (numTargets(thisTrial) == 2 || targetSide(thisTrial) == 1)
                % Draw cue on left
                Screen('FrameOval', ptbWindow, cueVal, theseCueRects(:,1), cueWidth_pix, cueWidth_pix);
            end

            if (thisItem == thisTarget(1)) && (numTargets(thisTrial) == 2 || targetSide(thisTrial) == 2)
                % Draw cue on right
                Screen('FrameOval', ptbWindow, cueVal, theseCueRects(:,2), cueWidth_pix, cueWidth_pix);
            end

            % Mark drawing ops as finished
            Screen('DrawingFinished', ptbWindow);

            % Done. Draw item.
            Screen('Flip', ptbWindow, startTime + (thisItem-1)*itemBlankDuration - flipProportion*flipInterval);

            % Blank screen with fixation

            if fixationOn&&(thisItem<nItems)
                Screen('FillOval', ptbWindow, fixVal, fixRect);
            end

            Screen('Flip', ptbWindow, startTime + (thisItem-1)*itemBlankDuration + itemDuration - flipProportion*flipInterval);
        end

        % Restore priority
        Priority(oldPriority);
                
        WaitSecs(stimResponseWait);

        % Draw response screen
        Screen('TextSize', ptbWindow, responsePoints);
        SetMouse(screenCentreX,screenCentreY, screenID);
        ShowCursor('Arrow');

        % Randomly determine response order
        if targetSide(thisTrial)==0
            responseOrder = randi(2); % 1=L first, 2=R first
        elseif targetSide(thisTrial)==1
            responseOrder = 1;
        else 
            responseOrder = 2;  
        end
        
        responseTime = NaN(1,2);
        selectedOrientation = NaN(1,2);
        selectedLetters = zeros(nLetters,2);

        rTimeStart = 0;
        
        for thisResponse = 1:2
            
            if ~(numTargets(thisTrial)==1 && thisResponse==2) % Don't loop twice if single target
            
                if ((thisResponse==1)&&(responseOrder==1))||((thisResponse==2)&&(responseOrder==2))
                    % Left/Top
                    responseSide = 1;
                    colValL = whiteVal;
                    colValR = dimVal;
                    validRects1 = validRespRects(1,:);
                    validRects2 = validRespRects(2,:);        
                else
                    % Right/Bottom
                    responseSide = 2;
                    colValL = dimVal;
                    colValR = whiteVal;
                    validRects1 = validRespRects(3,:); 
                    validRects2 = validRespRects(4,:);
                 end

                goodResponse = 0;

                while goodResponse==0

                    % Check that mouse is still within ptbWindow
                    [xMPos,yMPos] = GetMouse(screenID);

                    if xMPos < 10
                        SetMouse(screenCentreX, yMPos, screenID);
                        ShowCursor('Arrow');
                    end

                    % Draw attribute matrix
                    for thisLetter = 1:nLetters
                        if numTargets(thisTrial)==2 || targetSide(thisTrial)==1
                            %Left line-up
                            DrawFormattedText(ptbWindow, letterArray(thisLetter), (thisXresponse(1,thisLetter)-1.5*responsePoints), thisYresponse(1,thisLetter)+responsePoints, colValL);
                            DrawFormattedText(ptbWindow, letterArray(thisLetter), thisXresponse(1,thisLetter), thisYresponse(1,thisLetter)-responsePoints, colValL, [], true, true);
                        end
                        if numTargets(thisTrial)==2 || targetSide(thisTrial)==2
                            %Right line-up
                            DrawFormattedText(ptbWindow, letterArray(thisLetter), thisXresponse(2,thisLetter), thisYresponse(2,thisLetter)+responsePoints, colValR);
                            DrawFormattedText(ptbWindow, letterArray(thisLetter), (thisXresponse(2,thisLetter)+1.5*responsePoints), thisYresponse(2,thisLetter)-responsePoints, colValR, [], true, true);
                        end
                    end

                    % Draw response letters
                    Screen('TextSize', ptbWindow, letterPoints);
                    Screen('TextFont', ptbWindow, letterFont);

                    if sum(selectedLetters(:,1))
                        %Screen('DrawText', ptbWindow, letterArray(find(selectedLetters(:,1))), letterX(1), letterY(1), colValL, greyVal); %#ok<FNDSB>
                        %flipLetter = letterOrientaion(find(letterStream(1,:)==find(selectedLetters(:,1))));
                        DrawFormattedText(ptbWindow, letterArray(find(selectedLetters(:,1))), letterX(1), 'center', colValL, [], selectedOrientation(1), selectedOrientation(1)); %#ok<FNDSB>
                    elseif targetSide(thisTrial)~=2
                         DrawFormattedText(ptbWindow, letterNull, letterX(1), 'center', colValL);
                    end

                    if sum(selectedLetters(:,2))
                        %Screen('DrawText', ptbWindow, letterArray(find(selectedLetters(:,2))), letterX(2), letterY(2), colValR, greyVal); %#ok<FNDSB>
                        %flipLetter = letterOrientaion(find(letterStream(2,:)==find(selectedLetters(:,2))));
                        DrawFormattedText(ptbWindow, letterArray(find(selectedLetters(:,2))), letterX(2), 'center', colValR, [], selectedOrientation(2), selectedOrientation(2)); %#ok<FNDSB> %#ok<FNDSB>
                     elseif targetSide(thisTrial)~=1
                         DrawFormattedText(ptbWindow, letterNull, letterX(2), 'center', colValR);
                    end

                    Screen('TextSize', ptbWindow, responsePoints);

                    if sum(selectedLetters(:,responseSide)) % If attribute is selected
                        Screen('DrawTexture', ptbWindow, colourTex, [], goRect, [], [], [], [], []);
                        Screen('TextSize', ptbWindow, okPoints);
                        DrawFormattedText(ptbWindow, 'OK', 'center', 'center', blackVal);
                        Screen('TextSize', ptbWindow, responsePoints);
                    else
                        % Possibly draw black GO texture
                        % Screen('DrawTexture', ptbWindow, colourTex, [], goRect, [], [], [], [blackVal blackVal blackVal], []);
                        % DrawFormattedText(ptbWindow, 'OK', 'center', 'center', whiteVal);
                    end

                    if ~ rTimeStart
                        rTimeStart = Screen('Flip', ptbWindow);
                    else
                        Screen('Flip', ptbWindow);
                    end

                    % Get a click
                    [nClicks,xClick,yClick,buttonClick] = GetClicks(ptbWindow,0);
                    clickTime = GetSecs;

                    xClick = screenWidth-xClick;

                    % Check whether the click is in a letter rect  
                    if IsInRect(xClick,yClick,validRects1)||IsInRect(xClick,yClick,validRects2)         
                        [minDist,clickedLetter] = min(abs(yResponseArrayLetters_Centre-yClick));
                        selectedLetters(:,responseSide) = zeros(nLetters,1);
                        selectedLetters(clickedLetter,responseSide) = 1;
                        if IsInRect(xClick,yClick,validRects1)
                                selectedOrientation(responseSide) = 0;
                        else
                                selectedOrientation(responseSide) = 1;
                        end
                    end

                    if IsInRect(xClick,yClick,goRectResp)
                        if sum(selectedLetters(:,responseSide)) % If all attributes are selected
                            responseTime(responseSide) = clickTime-rTimeStart;
                            goodResponse=1;
                            Screen('Flip', ptbWindow);
                        else
                            % Possible error tone
                        end
                    end
                end
            end
        end
       
        % Store data for this trial
        allLetterOrder(thisTrial,:,:) = letterStream;
        allLetterOrientations(thisTrial,:) = letterOrientaion;
        
        
        if targetSide(thisTrial) == 1
            targetSaveSide = 1;
        elseif targetSide(thisTrial) == 2
            targetSaveSide = 2;
        else
            targetSaveSide = 1:2;
        end
            
        allTargets(thisTrial,targetSaveSide) = thisTarget(targetSaveSide);
        allTargetOrientations(thisTrial,targetSaveSide) = letterOrientaion(thisTarget(targetSaveSide));
        
        if targetSide(thisTrial) == 0
            allResponses(thisTrial,:) = [find(selectedLetters(:,1)) find(selectedLetters(:,2))];
        else
            allResponses(thisTrial,targetSaveSide) = find(selectedLetters(:,targetSaveSide));
        end
        
        allResponseOrientations(thisTrial,:) = selectedOrientation;
        allResponseOrder(thisTrial) = responseOrder;
        
        allRTs(thisTrial,:) = responseTime;
        
        if isempty(find(stopTrials==thisTrial, 1))&&(thisTrial~=totalTrials)
            Screen('TextSize', ptbWindow, instructionPoints);
            Screen('TextFont', ptbWindow, instructionFont);
            nextTrialInfoText = ['Click the mouse button to start the next trial.\n\n(Trial ' num2str(thisTrial+1) ' of ' num2str(totalTrials) ')'];
            DrawFormattedText(ptbWindow, nextTrialInfoText, 'center', 'center', whiteVal);
            HideCursor;
            Screen('Flip', ptbWindow);

            % Wait for click to continue
            GetClicks(ptbWindow,0);
            
        elseif (find(stopTrials==thisTrial, 1))~=nBlocks
            
            % Break
            thisBlock = find(stopTrials==thisTrial);
            remainingBlocks = nBlocks-thisBlock;
                
            if remainingBlocks == 1
                nextBlockInfoText = ['Thank you. The block is complete.'...
                        '\n\nPlease take a rest, then click the mouse button'...
                        '\n\nwhen you are ready to continue the experiment.\n\n'...
                        '(1 block remaining)'];
            else
                nextBlockInfoText = ['Thank you. The block is complete.'...
                        '\n\nPlease take a rest, then click the mouse button'...
                        '\n\nwhen you are ready to continue the experiment.\n\n'...
                        '\n\n(' num2str(remainingBlocks) ' blocks remaining)'];
            end

            Screen('TextSize', ptbWindow, instructionPoints);
            Screen('TextFont', ptbWindow, instructionFont);
            DrawFormattedText(ptbWindow, nextBlockInfoText, 'center', 'center', whiteVal);
            
            HideCursor;

            Screen('Flip', ptbWindow);

            WaitSecs(2);
            GetClicks([],0);
            Screen('Flip', ptbWindow);

            Screen('TextSize', ptbWindow, letterPoints);
            Screen('TextFont', ptbWindow, letterFont);
        end
        
    end

    
% Save all data
endTime = now;
oldDirectory = cd(saveDirectory);

% Make new filename
fileNumber = 1;
newFileName = [participantID '_' datestr(now,'YY-mm-DD') '_' num2str(fileNumber) '.mat'];

while exist([saveDirectory newFileName],'file')
    fileNumber = fileNumber+1;
    newFileName = [participantID '_' datestr(now,'YY-mm-DD') '_' num2str(fileNumber) '.mat'];
end

save(newFileName, 'participantID', 'randomSeed', 'itemRate', 'letterSeparation', 'letterEccentricity', 'allConditions', 'numTargets', 'targetSide', 'allLetterOrder', 'allLetterOrientations', 'allTargets', 'allTargetOrientations', 'allRTs', 'allResponses', 'allResponseOrientations', 'allResponseOrder', 'endTime');

% Update user data
cd(userDirectory);
nextSession = nextSession+1;
totalTime = [totalTime GetSecs-startTotalTime]; 
save(userDataFile,'sessionOrder','nextSession','totalTime');
cd(oldDirectory);

% Close onscreen window and audio port, release all ressources
WaitSecs(1);
Screen('TextSize', ptbWindow, instructionPoints);
Screen('TextFont', ptbWindow, instructionFont);
DrawFormattedText(ptbWindow, endText, 'center', 'center', whiteVal, textWrap);
Screen('Flip', ptbWindow);
WaitSecs(2);
GetClicks([],0);
Screen('Flip', ptbWindow);
PsychPortAudio('Close', ptbAudioPort);
Screen('CloseAll');