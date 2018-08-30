% Dual RSVP stream demo code
% PTG 11/07/12
Screen('CloseAll');
clear all; 

addpath(genpath('/Users/supalab/Documents/MATLAB/BackwardsLetters/E2/'));

AssertOpenGL;

% Initialise experiment parameters
saveDirectory = '/Users/supalab/Documents/MATLAB/BackwardsLetters/E2/Data/';
userDirectory = '/Users/supalab/Documents/MATLAB/BackwardsLetters/E2/UserData/';
nTrials = 100;            % 100, 50 (per condition)
nPracticeTrials = 16;    % 15, 10
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


% Get practise order

while true
    praticeOrder = input('Enter first practise block condition 1 or 2 (1 = Canonical, 2 = Inverted): ');
    if any(praticeOrder == 1:2)
        break
    end
end


% Update practice order and intructions

if praticeOrder == 1
    firstPracticeTxt =  ' in their normal upright orientation.\n\n';
    nextPracticeTxt =  ' rotated 180 degrees.\n\n';
    
else
    firstPracticeTxt = ' rotated 180 degrees.\n\n';  
    nextPracticeTxt =  ' in their normal upright orientation.\n\n';
    
end


oldDirectory = cd(userDirectory);
userDataFile = [participantID '_RotatedDualRSVP_Exp1_UserData.mat'];

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
letterArray = [char(65:71) char(74:76) char(81:82) char(84:86) char(89)]; %ABCDEFGJKLPQRTUVY %[char(65:66) char(68:85) char(87:90)];
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
expDegrees = [0 180];
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
instructionText1 = ['On each trial, a rapid stream of letters ' ...
    'will appear on the left and right of the screen.\n\n' ...
    'The two streams will be simultaneous, but the order of ' ...
    'letters will be different in each stream.\n\n' ...
    'Within each stream some of the letters will be rotated 180 degrees.\n\n' ...
    'Click the mouse button for more instructions.'];

instructionText2 = ['At some time during the trial, white circles will appear ' ...
    'around both streams.\n\n' ...
    'Your task is to report the letter and orientation that was ' ...
    'displayed in each stream at the time that the circles appeared.\n\n' ...
    'The letters may be any from a set of 16 letters from the alphabet.\n\n' ...
    'Click the mouse button for more instructions.'];

instructionText3 = ['Please keep your eyes fixated on the white dot at ' ...
    'the centre of the screen at all times during the trial.\n\n' ...
    'After each trial, please select your response by clicking on the letter and orientation ' ...
    'that you saw in each stream when the circles appeared.\n\n' ...
    'If you are unsure, please make your best guess.\n\n' ...
    'Click the mouse button for more instructions.'];

instructionText4 = ['If you have any questions, please ask the experimenter now.\n\n' ...
    'In the first practice block, the circled letters will always be' firstPracticeTxt ...
    'Click the mouse button to begin the practice trials.'];
    
endText = ['The experiment is complete.\n\n' ...
    'Please click the mouse button to end the program.'];

nextBlockConditionText = {'In the next block the oriention of circled letters will be random.\n\n', ...
   'In the next block the orientation of circled letters will be random.\n\n'}; 

firstBlockConditionText = {'For the remainder of the experiment \norientation of circled letters will be random.\n\n'};

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
%validRespRects(2,:) = CenterRectOnPoint([0 0 hSubtenseResponse_pix+responseLetterRect(3) responseLetterRect(4)], screenCentreX, yResponseArray(9));
%validRespRects(4,:) = CenterRectOnPoint([0 0 hSubtenseResponse_pix+responseLetterRect(3) responseLetterRect(4)], screenCentreX, yResponseArray(1)); %Right-outer (inveted)

xResponseArray = xResponseArray-[responseLetterRect(3) 0 0 letterRect(3) 0 letterRect(3) 0 0 responseLetterRect(3)]./2;
yResponseArray = yResponseArray-[responseLetterRect(4) 0 0 letterRect(4) 0 letterRect(4) 0 0 responseLetterRect(4)]./2;

goRect = CenterRectOnPoint(colourRect, screenCentreX, screenCentreY);

respShrinkPix = [shrinkFactor_pix shrinkFactor_pix -shrinkFactor_pix -shrinkFactor_pix];
goRectResp = goRect+respShrinkPix;

% Condition order

% Randomise trial/condition order
conditionList = repmat(1:2,1,nTrials); % Generate list of trial conditions, equal num per condition
allConditions = conditionList(randperm(length(conditionList))); % Make trial-condition random

% --------------------------  
    % Calculate temporal parameters
    itemRate = itemRates(1);                 % Item presentation rate (items per second)
    itemBlankDuration = 1.0/itemRate;                                       % Duration in secs of item + blank
    itemDuration = (itemBlankRatio/(itemBlankRatio+1))*itemBlankDuration;   % Duration in secs of item
    blankDuration = itemBlankDuration-itemDuration;    


% Instruction screen
Screen('TextSize', ptbWindow, instructionPoints);
Screen('TextFont', ptbWindow, instructionFont);

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

    DrawFormattedText(ptbWindow, instructionText4, 'center', 'center', whiteVal, textWrap);
    WaitSecs(2);
    GetClicks([],0);
    
    Screen('Flip', ptbWindow);

    WaitSecs(2);
    GetClicks([],0);
    Screen('Flip', ptbWindow);

    Screen('TextSize', ptbWindow, letterPoints);
    Screen('TextFont', ptbWindow, letterFont);
    
    startBlockTime = GetSecs;
    clickTime = startBlockTime;
    HideCursor;
    
    
    for thisPracticeBlock = 1:2
            
        for thisTrial = 1:nPracticeTrials

            Screen('TextSize', ptbWindow, letterPoints);
            Screen('TextFont', ptbWindow, letterFont);
            
            % Set practise block condtion
            if ((praticeOrder==1)&&(thisPracticeBlock==1))||((praticeOrder==2)&&(thisPracticeBlock==2))
                thisCondition = 1;
            else
                thisCondition = 2;
            end
            
            if simulTargets
                thisTarget = repmat(randi(nPossTargets)+nLeadTrailItems,1,2);
            else
                thisTarget = randi(nPossTargets,[1 2])+nLeadTrailItems;
            end

            letterStream = mod([randperm(nItems); randperm(nItems)],nLetters)+1;
            letterOrientaion =  repmat(0:1,1,nLetters/2); % Must be an even number of letters in letterArray
            %AH Now shuffle that list of 0's and 1's, by creating randomly shuffled indices
            letterOrientaion = letterOrientaion(randperm(length(letterOrientaion)));
            
            %Check target has correct orientation
            if (thisCondition==1)&&(letterOrientaion(thisTarget(1,1))==1)
            	%AH I presume the above being true means target has the wrong orientation
                uprightLetters = find(letterOrientaion==0); %AH uprightLetters now is a list of indices of where the stream has the correct orientation 
                %AH Now, shuffle the 
                %AH What does the ",1" do in randperm
                %AH Seems like its setting a random upright letter to upside down,
                %AH that way can then set the actual target to rightside up
                %AH randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n inclusive.  In other words, picks one of the upright letter positions
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
            %Screen('DrawText', ptbWindow, letterArray(letterStream(1,1)), letterX(1), letterY(1), letterVal, greyVal);
            %Screen('DrawText', ptbWindow, letterArray(letterStream(2,1)), letterX(2), letterY(2), letterVal, greyVal); 
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
                % Screen('DrawText', ptbWindow, letterArray(letterStream(1,thisItem)), letterX(1), letterY(1), letterVal, greyVal);
                % Screen('DrawText', ptbWindow, letterArray(letterStream(2,thisItem)), letterX(2), letterY(2), letterVal, greyVal);
                DrawFormattedText(ptbWindow, letterArray(letterStream(1,thisItem)), letterX(1), 'center', letterVal, [], letterOrientaion(thisItem), letterOrientaion(thisItem));
                DrawFormattedText(ptbWindow, letterArray(letterStream(2,thisItem)), letterX(2), 'center', letterVal, [], letterOrientaion(thisItem), letterOrientaion(thisItem));
                
                if fixationOn
                    Screen('FillOval', ptbWindow, fixVal, fixRect);
                end

                if thisItem == thisTarget(1)
                    % Draw cue on left
                    Screen('FrameOval', ptbWindow, cueVal, theseCueRects(:,1), cueWidth_pix, cueWidth_pix);
                end

                if thisItem == thisTarget(2)
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
            responseOrder = randi(2); % 1=L first, 2=R first
            responseTime = NaN(1,2);
            selectedLetters = zeros(nLetters,2);
            selectedOrientation = NaN(1,2);
            rTimeStart = 0;

            for thisResponse = 1:2
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
                    
                      % Match oreintation to that in stream
%                     if letterOrientaion(find(letterStream(1,:)==thisLetter))==0
%                         adjustYPos1 = responsePoints;
%                         flipLetter1 = false;
%                     else
%                         adjustYPos1 = -responsePoints;
%                         flipLetter1 = true;
%                     end
%                     
%                     if letterOrientaion(find(letterStream(2,:)==thisLetter))==0
%                         adjustYPos2 = responsePoints;
%                         flipLetter2 = false;
%                     else
%                         adjustYPos2 = -responsePoints;
%                         flipLetter2 = true;
%                     end
                              
                    %Screen('DrawText', ptbWindow, letterArray(thisLetter), thisXresponse(1,thisLetter), thisYresponse(1,thisLetter), colValL, greyVal);
                    %Screen('DrawText', ptbWindow, letterArray(thisLetter), thisXresponse(2,thisLetter), thisYresponse(2,thisLetter), colValR, greyVal);
                    
                    %Left line-up
                    DrawFormattedText(ptbWindow, letterArray(thisLetter), (thisXresponse(1,thisLetter)-1.5*responsePoints), thisYresponse(1,thisLetter)+responsePoints, colValL);
                    DrawFormattedText(ptbWindow, letterArray(thisLetter), thisXresponse(1,thisLetter), thisYresponse(1,thisLetter)-responsePoints, colValL, [], true, true);
                    %Right line-up
                    DrawFormattedText(ptbWindow, letterArray(thisLetter), thisXresponse(2,thisLetter), thisYresponse(2,thisLetter)+responsePoints, colValR);
                    DrawFormattedText(ptbWindow, letterArray(thisLetter), (thisXresponse(2,thisLetter)+1.5*responsePoints), thisYresponse(2,thisLetter)-responsePoints, colValR, [], true, true);
                
                end

                    % Draw response letters
                    Screen('TextSize', ptbWindow, letterPoints);
                    Screen('TextFont', ptbWindow, letterFont);

                    if sum(selectedLetters(:,1))
                        %Screen('DrawText', ptbWindow, letterArray(find(selectedLetters(:,1))), letterX(1), letterY(1), colValL, greyVal); %#ok<FNDSB>
                        %flipLetter = letterOrientaion(find(letterStream(1,:)==find(selectedLetters(:,1)))
                        DrawFormattedText(ptbWindow, letterArray(find(selectedLetters(:,1))), letterX(1), 'center', colValL, [], selectedOrientation(1), selectedOrientation(1)); %#ok<FNDSB>
                    else
                        DrawFormattedText(ptbWindow, letterNull, letterX(1), 'center', colValL);
                    end

                    if sum(selectedLetters(:,2))
                        %Screen('DrawText', ptbWindow, letterArray(find(selectedLetters(:,2))), letterX(2), letterY(2), colValR, greyVal); %#ok<FNDSB>
                        %flipLetter = letterOrientaion(find(letterStream(2,:)==find(selectedLetters(:,2))));
                        DrawFormattedText(ptbWindow, letterArray(find(selectedLetters(:,2))), letterX(2), 'center', colValR, [], selectedOrientation(2), selectedOrientation(2)); %#ok<FNDSB>
                    else
                        DrawFormattedText(ptbWindow, letterNull, letterX(2), 'center', colValL);
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


            if thisTrial < nPracticeTrials
                nextTrialInfoText = ['Click the mouse button to start the next practice trial.\n\n(Trial ' num2str(thisTrial+1) ' of ' num2str(nPracticeTrials) ')'];
                Screen('TextSize', ptbWindow, instructionPoints);
                Screen('TextFont', ptbWindow, instructionFont);
                DrawFormattedText(ptbWindow, nextTrialInfoText, 'center', 'center', whiteVal);
                HideCursor;
                Screen('Flip', ptbWindow);

                % Wait for click to continue
                GetClicks(ptbWindow,0);

            elseif thisPracticeBlock==1                
                
                nextTrialInfoText = ['In the next block of practice trials, the circled letters\n will always be' nextPracticeTxt ...                           
                'Click the mouse button to start the next practice trial.\n\n' ...
                '(Trial ' num2str(1) ' of ' num2str(nPracticeTrials) ')'];
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
                firstBlockConditionText{allConditions(1)} ...
                'Please click the mouse button to begin the experiment.'];
                
                Screen('TextSize', ptbWindow, instructionPoints);
                Screen('TextFont', ptbWindow, instructionFont);
                DrawFormattedText(ptbWindow, nextTrialInfoText, 'center', 'center', whiteVal);
            end

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
        %Screen('DrawText', ptbWindow, letterArray(letterStream(1,1)), letterX(1), letterY(1), letterVal, greyVal);
        %Screen('DrawText', ptbWindow, letterArray(letterStream(2,1)), letterX(2), letterY(2), letterVal, greyVal);
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
            % Screen('DrawText', ptbWindow, letterArray(letterStream(1,thisItem)), letterX(1), letterY(1), letterVal, greyVal);
            % Screen('DrawText', ptbWindow, letterArray(letterStream(2,thisItem)), letterX(2), letterY(2), letterVal, greyVal);
            DrawFormattedText(ptbWindow, letterArray(letterStream(1,thisItem)), letterX(1), 'center', letterVal, [], letterOrientaion(thisItem), letterOrientaion(thisItem));
            DrawFormattedText(ptbWindow, letterArray(letterStream(2,thisItem)), letterX(2), 'center', letterVal, [], letterOrientaion(thisItem), letterOrientaion(thisItem));

            if fixationOn
                Screen('FillOval', ptbWindow, fixVal, fixRect);
            end

            if thisItem == thisTarget(1)
                % Draw cue on left
                Screen('FrameOval', ptbWindow, cueVal, theseCueRects(:,1), cueWidth_pix, cueWidth_pix);
            end

            if thisItem == thisTarget(2)
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
        responseOrder = randi(2); % 1=L first, 2=R first
        responseTime = NaN(1,2);
        selectedOrientation = NaN(1,2);
        selectedLetters = zeros(nLetters,2);

        rTimeStart = 0;
        
        for thisResponse = 1:2
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
                
                       % Match oreintation to that in stream
%                      if letterOrientaion(find(letterStream(1,:)==thisLetter))==0
%                         adjustYPos1 = responsePoints;
%                         flipLetter1 = false;
%                     else
%                         adjustYPos1 = -responsePoints;
%                         flipLetter1 = true;
%                     end
%                     
%                     if letterOrientaion(find(letterStream(2,:)==thisLetter))==0
%                         adjustYPos2 = responsePoints;
%                         flipLetter2 = false;
%                     else
%                         adjustYPos2 = -responsePoints;
%                         flipLetter2 = true;
%                     end
                    
                    %Screen('DrawText', ptbWindow, letterArray(thisLetter), thisXresponse(1,thisLetter), thisYresponse(1,thisLetter), colValL, greyVal);
                    %Screen('DrawText', ptbWindow, letterArray(thisLetter), thisXresponse(2,thisLetter), thisYresponse(2,thisLetter), colValR, greyVal);
                    
                    %Left line-up
                    DrawFormattedText(ptbWindow, letterArray(thisLetter), (thisXresponse(1,thisLetter)-1.5*responsePoints), thisYresponse(1,thisLetter)+responsePoints, colValL);
                    DrawFormattedText(ptbWindow, letterArray(thisLetter), thisXresponse(1,thisLetter), thisYresponse(1,thisLetter)-responsePoints, colValL, [], true, true);
                    %Right line-up
                    DrawFormattedText(ptbWindow, letterArray(thisLetter), thisXresponse(2,thisLetter), thisYresponse(2,thisLetter)+responsePoints, colValR);
                    DrawFormattedText(ptbWindow, letterArray(thisLetter), (thisXresponse(2,thisLetter)+1.5*responsePoints), thisYresponse(2,thisLetter)-responsePoints, colValR, [], true, true);
                    
                end
         
                % Draw response letters
                Screen('TextSize', ptbWindow, letterPoints);
                Screen('TextFont', ptbWindow, letterFont);

                if sum(selectedLetters(:,1))
                    %Screen('DrawText', ptbWindow, letterArray(find(selectedLetters(:,1))), letterX(1), letterY(1), colValL, greyVal); %#ok<FNDSB>
                    %flipLetter = letterOrientaion(find(letterStream(1,:)==find(selectedLetters(:,1))));
                    DrawFormattedText(ptbWindow, letterArray(find(selectedLetters(:,1))), letterX(1), 'center', colValL, [], selectedOrientation(1), selectedOrientation(1)); %#ok<FNDSB>
                else
                     DrawFormattedText(ptbWindow, letterNull, letterX(1), 'center', colValL);
                end

                if sum(selectedLetters(:,2))
                    %Screen('DrawText', ptbWindow, letterArray(find(selectedLetters(:,2))), letterX(2), letterY(2), colValR, greyVal); %#ok<FNDSB>
                    %flipLetter = letterOrientaion(find(letterStream(2,:)==find(selectedLetters(:,2))));
                    DrawFormattedText(ptbWindow, letterArray(find(selectedLetters(:,2))), letterX(2), 'center', colValR, [], selectedOrientation(2), selectedOrientation(2)); %#ok<FNDSB> %#ok<FNDSB>
                else
                     DrawFormattedText(ptbWindow, letterNull, letterX(2), 'center', colValL);
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
       
        % Store data for this trial
        allLetterOrder(thisTrial,:,:) = letterStream;
        allLetterOrientations(thisTrial,:) = letterOrientaion;
        
        allTargets(thisTrial,:) = thisTarget;
        allTargetOrientations(thisTrial,:) = letterOrientaion(thisTarget);
        
        allResponses(thisTrial,:) = [find(selectedLetters(:,1)) find(selectedLetters(:,2))];
        allResponseOrientations(thisTrial,:) = selectedOrientation;
        
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
                        nextBlockConditionText{allConditions(thisTrial+1)} ...
                        '(1 block remaining)'];
            else
                nextBlockInfoText = ['Thank you. The block is complete.'...
                        '\n\nPlease take a rest, then click the mouse button'...
                        '\n\nwhen you are ready to continue the experiment.\n\n'...
                        nextBlockConditionText{allConditions(thisTrial+1)} ...
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

save(newFileName, 'participantID', 'randomSeed', 'itemRate', 'letterSeparation', 'letterEccentricity', 'allLetterOrder', 'allLetterOrientations', 'allTargets', 'allTargetOrientations', 'allRTs', 'allResponses', 'allResponseOrientations', 'allConditions', 'endTime', 'expDegrees');

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