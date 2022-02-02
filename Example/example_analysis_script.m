%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE DYNAMIC ENDOCAST BY BONE AUTOMATION + RCVC CALCULATION
%
% This script walks you through how to use the dynamicEndocastByBone
% function in the command line (without GUIs), for automating across
% trials. It also demonstrates how to calculate and plot RCVC.
%
% The example data provided are from Whitlow et al. 2022, but feel free to
% adopt this script to your own data by replacing the filenames / bone
% names etc., to fit your own data.
%
% Written by J.D. Laurence-Chasen, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 0. Download example data
%
% Download the 4 example csv files from the "Example" folder of the
% repository.
%
% Note that these example files have inconsistent column/row headers--this
% is totally fine. The function with automatically detect whether there is
% a header row or frame number column with each individual file

%% Step 1. Set parameters and filenames

% Where are the data files, where do you want to save the output file?
% Replace with your own paths %
datapath = 'C:\Users\jdlc7\Documents\GitHub\dynamicEndocastByBone\Example\'; 
savepath = 'C:\Users\jdlc7\Documents\GitHub\dynamicEndocastByBone\Example\VolumeData';

refbone = 'Vol_NC'; % THIS MUST MATCH THE PREFIX IN THE LOCATOR FILE COLUMN HEADERS
freezeIncrement = 10; % See Whitlow et al. 2022--experiment with different increments
filter = true; 
objFolder = []; % DO NOT save hull OBJs

trials = {'P1_Trial01', 'P1_Trial02'}; % create a list of trials
% make sure all trial files (xyz locator, RBT) start with the same prefix

%% Step 2. Run dynamicEndocastByBone

endocast_data = {}; % Output variable 

for t = 1:length(trials) % for every trial
    
    % Get file names
    % (Assumes your locator files are in the format of:
    % 'P1_Trial01_VolumeLocators.csv')
    % Change to match your file naming conventions
    
    XYZfile = [trials{t} '_VolumeLocators.csv']; % Locators from Maya
    RBTfile = [trials{t} '_NeurocraniumRBT.csv']; % Reference RBT from XMALab
    savefile = [trials{t} '_VolumeData.csv']; % Name of output file (you pick)
 
    % Do it
    endocast_data{t} = ...
        dynamicEndocastByBone(datapath, XYZfile, RBTfile, refbone, ...
        freezeIncrement,filter,savepath,savefile,objFolder);
end

%% Step 3. Calculate RCVC

rcvc = {};
bones = {'CH' 'LJ' 'OP' 'SP'}; % Bones to calculate RCVC for

for t = 1:length(trials)
    
    % Which frames have data (note that there are NaNs in output file)
    frames = ~isnan(endocast_data{t}.DeltaVolumeFull);
    
    impacts = zeros(sum(frames),1); % pre-allocate impacts variable
    
    % Calculate impacts (Delta Full - Delta Frozen Bone)
    for b = 1:length(bones) 
        impacts(:,b) = endocast_data{t}.DeltaVolumeFull(frames) - ...
            endocast_data{t}.(['DeltaVolumeFrozen_Vol_' bones{b}])(frames);
    end
    
    sum_of_impacts = sum(abs(impacts),2);
    
    for b = 1:length(bones)
        
        % Run equation
        rcvc{t}(:,b) = impacts(:,b) ./ sum_of_impacts;
        
    end
end

%% Step 4. Plot

figure;
for b = 1:length(bones)
    subplot(1,length(bones),b)
    hold on
    for t = 1:length(trials)
        plot(rcvc{t}(:,b))
    end
    legend(trials,'Interpreter', 'none')
    title([bones(b) ' RCVC'])
    xlabel('Frame')
    ylabel('RCVC')
end





