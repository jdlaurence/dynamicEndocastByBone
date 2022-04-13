%% dynamicEndocastByBone

% This script is built off of dynamicEndocast (by Ariel Camp and Peter
% Falkingham), and calculates the specific contribution of indvidual bones
% to the change in volume of the endocast. It works by 'freezing' the
% motion of individual bones (or regions of locators) relative to a
% reference bone (usually the cranium).
%
% The output of the script is a file that has the derivative of volume with
% all locators, followed by the derivative of volume with specific bones
% frozen for a given frame (or running frame range, if averaging is used).
% The difference between the full volume and frozen volume at a time point
% is the 'contribution' of that bone's motion to the overall volume change.
%
% Importantly, it identifies the bones/regions by the columm headers in the
% locator XYZ file. So if you want locators to be grouped in a region, be
% sure their names follow the format: locatornameXX (must end in a number).
%
% Inputs
%   datapath: string - where are files located 
%   XYZfile: string - path of XYZ locator coordinates (output of maya)
%   RBTfile: string - path of reference bone RBT file (from XMALab)
%   refbone: string - name of reference bone in Column names of XYZ file
%   freezeIncrement: double - >0 integer of freeze increment
%   savepath: string - where new files should be saved
%   savefile: string - filename of output 
%   objFolderPath: string (optional - where to save obj files
%   bridges: cell array: Names of locators to bridge (set to [] for none))
%
% Example usage:
%
%   With GUIs
%       dynamicEndocastByBone;
%   Without GUIs, no OBJ save
%       dynamicEndocastByBone(datapath,XYZfile,RBTfile,refbone,freezeIncrement,savepath,savefile,[]);
%   Without GUIs, with OBJ save
%       dynamicEndocastByBone(datapath,XYZfile,RBTfile,refbone,freezeIncrement,savepath,savefile,objFolderPath);
%
% Written by J.D. Laurence-Chasen 8/5/2020
%
% History: 
% Added command line support/made function and fixed frame averaging
% 11/23/2020
%
% Add volume export functionality
% 2/18/2021

function [dataout] = dynamicEndocastByBone(datapath,XYZfile,RBTfile,refbone,freezeIncrement,filter,savepath,savefile,objFolder,bridges,varargin)

%% Parse Inputs And Run GUIs if necessary

if nargin == 0 % GUIs
    
    % % 1.  Read Locator file % % 
    [FileName,PathName] = uigetfile('*.csv','Select volume xyz coordinates csv file');
    cd(PathName);
    
    % % 2. Read Reference RBT % % 
    %             (assumes default output settings from XMALab)
    [RBTfile] = uigetfile('*.csv','Select reference bone rigid body transformation file');
    try % no header
        rbt = csvread(RBTfile,0,0);
    catch % header
       rbt = csvread(RBTfile,1,0); 
    end
    if size(rbt,2) == 17 % there are frame numbers
        rbt = rbt(:,2:end);
        multRBTs = false;
    elseif size(rbt,2) == 16
        multRBTs = false;
    elseif size(rbt,2) > 17 % There are multiple RBTs in the ref rbt
        if rem(size(rbt,2),16) == 1
            % there's a frames column and multiple RBTs
            rbt = rbt(:,2:end);
            multRBTs = true;
        elseif rem(size(rbt,2),16) == 0
            multRBTs = true;
        else
            error('Please use default RBT export settings from XMALab (no frames column)')
        end
    end
    
    % % 3. Try to automatically identify bones from row headers % % 
    % I'm going to assume the names are alpha characters followed by a
    % number followed by _tx/_ty/_tz
    XYZfile = readtable(FileName,'delimiter','comma');
    varnames = XYZfile.Properties.VariableNames; % get column names
    XYZfile = csvread(FileName,1,0); % re-read file as array (instead of table)
    if XYZfile(:,end) == 0 %check if there's an extra column of zeros
        XYZfile = XYZfile(:,1:end-1);
    end
    % Are there frame numbers?
    hasFN_XYZ = CheckForFrameNumbers(XYZfile);
    if hasFN_XYZ
        varnames = varnames(2:end)'; % first var is Frame
    elseif ~hasFN_XYZ
        varnames = varnames';
    else
        error('Irregular number of columns in XYZ point file')
    end
    croppedvars = varnames;
    for v = 1:length(varnames)
        croppedvars{v,1} = varnames{v}(1:end-3); % crop off _tx/_ty/_tz
    end
    varnamesnew = unique(croppedvars,'stable'); % get rid of duplicates
    % find numbers at end of cell
    nums = regexp(varnamesnew,'\d*','Match');
    varnames_nonums = cell(length(nums),1);
    for v = 1:length(nums)
        if ~isempty(nums{v})
            tmpv = varnamesnew{v}(1:end-length(nums{v}{1}));
            varnames_nonums{v,1} = tmpv;
        end
    end
    
    % % 4. Select reference bone from identified bones % % 
    bones = unique(varnames_nonums(~cellfun(@isempty,varnames_nonums)),'stable');
    [indx,tf] = listdlg('PromptString',{'Select the reference bone corresponding to the RBT file you just selected',''...
        'Make sure this list looks correct.',''},...
        'SelectionMode','single','ListString',bones);
    switch tf
        case 1
            refbone = bones{indx};
        otherwise
            error('you gotta pick a reference bone!!!')
    end
    
    bones(indx) = []; % get rid of reference bone
    
    % % 2.2 Isolate reference RBT, if there are multiple RBTs
    if multRBTs      
        cols = 16*(1-indx) + 1: 16*(1-indx)+16;
        rbt = rbt(:,cols);      
    end
    
    
    % % 5. Get column indicies of locators grouped by bone/region % % 
    reflocpts = find(contains(varnamesnew,refbone));
    refloccols = find(contains(varnames,refbone));
    bonelocpts = {};
    boneloccols = {};
    for b = 1:length(bones)
        bonelocpts{b} = find(contains(varnamesnew,bones{b}));
        boneloccols{b} = find(contains(varnames,bones{b}));
    end
    
    % % 6. Get freeze increment info % % 
    prompt = {'How many frames do you want the freeze increment to be?'};
    dlgtitle = 'Freeze Increment';
    definput = {'10'};
    freezeIncrement = inputdlg(prompt,dlgtitle,[1 50],definput);
    freezeIncrement = str2double(freezeIncrement{1});
    
    % % 7. Filter data? 
        filteranswer = questdlg('Filter data with a 3-frame moving average?', ...
	'Moving Average Filter', ...
	'Yes','No','Yes');
    % Handle response
    switch filteranswer
        case 'Yes'
            filter = true;
        case 'No'
            filter = false;
        otherwise
            filter = false;
            disp('Defaulting to NO filter')
    end
    
    % % 8. Get save path % % 
    [DataFile,DataPath] = uiputfile('*.csv','Save data');
    
    % % 9. Save obj files?  % % 
    answer = questdlg('Export volumes as .obj files?', ...
	'OBJ Export', ...
	'Yes','No','No');
% Handle response
switch answer
    case 'Yes'
        objexp = 1;
        
        % Get folder for save
        objFolder = uigetdir(PathName,'Select Folder for obj files');
        % Make save folders
        objfoldernames = {};
        objfoldernames{1} = fullfile(objFolder,'Reference');
        mkdir(objfoldernames{1});
        for b = 1:length(bones)
            objfoldernames{b+1} = fullfile(objFolder,['Frozen_' bones{b}]);
            mkdir(objfoldernames{b+1});
        end
        
        
    case 'No'
        objexp = 0;
    otherwise
        objexp = 0;
        disp('Defaulting to NO export')
end
    
elseif nargin ~= 10
    
    error('Must provide 9 input arguments (see example at start of function)')
    
else
    
    % % Read RBT
    PathName = datapath;
    cd(PathName);
    try % no header
        rbt = csvread(RBTfile,0,0);
    catch % header
        rbt = csvread(RBTfile,1,0);
    end
    if size(rbt,2) == 17 % there are frame numbers
        rbt = rbt(:,2:end);
        multRBTs = false;
    elseif size(rbt,2) == 16
        multRBTs = false;
    elseif size(rbt,2) > 17 % There are multiple RBTs in the ref rbt
        if rem(size(rbt,2),16) == 1
            % there's a frames column and multiple RBTs
            rbt = rbt(:,2:end);
            multRBTs = true;
        elseif rem(size(rbt,2),16) == 0
            multRBTs = true;
        else
            error('Please use default RBT export settings from XMALab (no frames column)')
        end
    end
    
    FileName = XYZfile;
    % Read Locator File
    XYZfile = readtable(FileName,'delimiter','comma');
    varnames = XYZfile.Properties.VariableNames; % get column names
    XYZfile = csvread(FileName,1,0); % re-read file as array (instead of table)
    if XYZfile(:,end) == 0 %check if there's an extra column of zeros
        XYZfile = XYZfile(:,1:end-1);
    end
    % Are there frame numbers?
    hasFN_XYZ = CheckForFrameNumbers(XYZfile);
    if hasFN_XYZ
        varnames = varnames(2:end)'; % first var is Frame
    elseif ~hasFN_XYZ
        varnames = varnames';
    else
        error('Irregular number of columns in XYZ point file')
    end
    % try to automatically identify bones
    % I'm going to assume the names are alpha characters followed by a
    % number followed by _tx/_ty/_tz
    
    croppedvars = varnames;
    for v = 1:length(varnames)
        croppedvars{v,1} = varnames{v}(1:end-3); % crop off _tx/_ty/_tz
    end
    varnamesnew = unique(croppedvars,'stable'); % get rid of duplicates
    
    % find numbers at end of cell
    nums = regexp(varnamesnew,'\d*','Match');
    varnames_nonums = cell(length(nums),1);
    for v = 1:length(nums)
        if ~isempty(nums{v})
            tmpv = varnamesnew{v}(1:end-length(nums{v}{1}));
            varnames_nonums{v,1} = tmpv;
        end
    end
    
    bones = unique(varnames_nonums(~cellfun(@isempty,varnames_nonums)),'stable');
    indx = find(strcmp(bones,refbone));
    if length(indx) ~= 1
        error('Make sure refbone corresponds to a bone/region name in the locator file')
    end
    
    bones(indx) = []; % get rid of reference bone
    
    if multRBTs      
        cols = 16*(1-indx) + 1: 16*(1-indx)+16;
        rbt = rbt(:,cols);      
    end
    
    
    % % Get column indicies of locators grouped by bone/region
    reflocpts = find(contains(varnamesnew,refbone));
    refloccols = find(contains(varnames,refbone));
    bonelocpts = {};
    boneloccols = {};
    for b = 1:length(bones)
        bonelocpts{b} = find(contains(varnamesnew,bones{b}));
        boneloccols{b} = find(contains(varnames,bones{b}));
    end
    
    % % Get locator bridge info % % 
    nbridges = size(bridges,1);
    % bridging settings
    dist_thresh = 0; % Minimum distance between generated locators
    max_pts = 10; % Maximum number of  generated locators
    
    % % Save info
    DataPath = savepath;
    DataFile = savefile;   
    
    if ~isempty(objFolder) % if saveing objs, make folders
        objexp = 1;
        objfoldernames = {};
        objfoldernames{1} = fullfile(objFolder,'Reference');
        if ~isfolder(objfoldernames{1})
            mkdir(objfoldernames{1});
        end
        for b = 1:length(bones)
            objfoldernames{b+1} = fullfile(objFolder,['Frozen_' bones{b}]);
            if ~isfolder(objfoldernames{b+1})
                mkdir(objfoldernames{b+1});
            end
        end
    else
       objexp = false; 
    end
    
end

%% Volumetric calculations

disp('Running endocast calculations...')

% Are there frame numbers?
if hasFN_XYZ %yes
    FrameNumbers = XYZfile(:,1);
    XYZpoints = XYZfile(:,2:end);
elseif ~hasFN_XYZ % no
    FrameNumbers = 1:size(XYZfile,1);
    FrameNumbers = FrameNumbers';
    XYZpoints = XYZfile;
end

%Number of frames:
nbones = length(bones);
nframes = size(XYZpoints,1);
npoints = size(XYZpoints,2)/3;
hullpoints = zeros(npoints,3);
alphavolume = zeros(nframes,1);
FinalVolumes = NaN(nframes,nbones*2+3);

for x = 1:nframes
    
    % Regular volume calculation %
    for y=1:npoints
        hullpoints(y,:) = XYZpoints(x,y*3-2:y*3);
    end
    [uniquepoints,~,~] = unique(hullpoints,'rows','stable'); %deletes duplicate XYZpoints
    
    % GENERATE BRIDGES (new option)
    if nbridges > 0
        for br = 1:nbridges
            % get indices for locator endpoints
            p1idx = find(strcmp(bridges{br,1},varnamesnew));
            p2idx = find(strcmp(bridges{br,2},varnamesnew));
            
            p1 = hullpoints(p1idx,:);
            p2 = hullpoints(p2idx,:);
            
            % calculate bridge
            midpoints = calcMidpoints(p1, p2, dist_thresh, max_pts);
            
            % add to locator set
            uniquepoints = [uniquepoints; midpoints];
                        
        end
    end
    
    n=2; %This is the value that determines the maximum radius of curvature
    shp = alphaShape(uniquepoints,n);
    totalvol = volume(shp);
    
    alphavolume(x,1) = totalvol;
    FinalVolumes(x,2) = totalvol;
    
    if objexp % save obj
        
        fileout = fullfile(objfoldernames{1}, ...
            ['Alpha_' num2str(FrameNumbers(x)) '_Hull.obj']);
        exportAlphaShapeAsOBJ(shp,FrameNumbers(x),fileout);
        
    end
    
    % The ~exciting~ bit %
    %LOGIC: In every frame, 'freeze' each of the bones sequentially,
    %relative to the reference bone, so they are in the pose from the
    %previous frame (technically, the pose from the nframeAvg previously.
    %Then recalculate the volume and see the difference between the new
    %volumne and original. This difference represents the contribution of
    %the motion of the bone in question to the overall volume change across
    %that nFrameAvg time period.
    
    
    if x > freezeIncrement % start freezing after freezeIncrement frames
        
        % indices
        curr_frame = x;
        prev_frame = x-freezeIncrement;
        store_frame = x-(floor(freezeIncrement/2));
        
        deltaFull = totalvol - alphavolume(prev_frame); % delta of unfrozen for increment
        FinalVolumes(store_frame,nbones+3) = deltaFull; % store in final array
        
        for b = 1:nbones
            
            XYZpointsFrozen = XYZpoints;
            
            pts = XYZpoints(prev_frame,boneloccols{b}); % points to freeze
            
            currentTm = reshape(rbt(x,:),[4 4]); % reference RBT for current frame
            prevTm = reshape(rbt(prev_frame,:),[4 4]); % reference RBT for previous frame
            newpts = [];
            
            for i = 1:length(pts)/3 % for each point in bone of interest,, 'freeze',
                %i.e. hold in place relative to reference bone by applying ref RBTs
                
                newPt = applyTm(pts((i-1)*3+1:(i-1)*3+3),inv(prevTm)); % put back in CT space
                newPt = applyTm(newPt,currentTm); % back into cube space, relative to ref
                newpts((i-1)*3+1:(i-1)*3+3) = newPt;
                
            end
            
            % replace bone points, with new, 'frozen' points
            XYZpointsFrozen(x,boneloccols{b}) = newpts;
            
            % recalculate alpha shape with new points
            for y=1:npoints
                hullpoints(y,:) = XYZpointsFrozen(x,y*3-2:y*3);
            end
            
            [uniquepoints,~,~] = unique(hullpoints,'rows','stable'); %deletes duplicate XYZpoints
            
            % GENERATE BRIDGES (new option)
            if nbridges > 0
                for br = 1:nbridges
                    % get indices for locator endpoints
                    p1idx = find(strcmp(bridges{br,1},varnamesnew));
                    p2idx = find(strcmp(bridges{br,2},varnamesnew));
                    
                    p1 = hullpoints(p1idx,:);
                    p2 = hullpoints(p2idx,:);
                    
                    % calculate bridge
                    midpoints = calcMidpoints(p1, p2, dist_thresh, max_pts);
                    
                    % add to locator set
                    uniquepoints = [uniquepoints; midpoints];
                    
                end
            end
            
            n=2; %This is the value that determines the maximum radius of curvature
            shp = alphaShape(uniquepoints,n);
            totalvolfrozen = volume(shp);
            FinalVolumes(store_frame,b+2) = totalvolfrozen;
            
            newdelta = totalvolfrozen - alphavolume(prev_frame);
            FinalVolumes(store_frame, nbones+3+b) = newdelta;
            
            if objexp % save obj
                
                  fileout = fullfile(objfoldernames{b+1}, ...
                      ['Alpha_' num2str(FrameNumbers(x)) '_Hull.obj']);
                  exportAlphaShapeAsOBJ(shp,FrameNumbers(x),fileout);          
            end            
        end
    end   
end


%% Filter if toggle on
if filter
    for i = 1:size(FinalVolumes,2)-1
        FinalVolumes(:,i+1) = movmean(FinalVolumes(:,i+1),3);
        
    end
end
        
%% file save
FinalVolumes(:,1) = FrameNumbers;

varNames = cell(1,nbones*2+3);
varNames{1} = 'FrameNumber';
varNames{2} = 'VolumeFull';
varNames{nbones+3} = 'DeltaVolumeFull';
for b = 1:nbones
    
    varNames{b+2} = ['VolumeFrozen_' bones{b}];
    varNames{nbones+3+b} = ['DeltaVolumeFrozen_' bones{b}];
    
end

dataout = array2table(FinalVolumes,'VariableNames',varNames);

% Create save path if it doesn't exist
if ~isfolder(DataPath)
    mkdir(DataPath)
end

cd(DataPath);
writetable(dataout,DataFile);
disp(['All done! Data saved to ' DataFile])

function [framenumbers] = CheckForFrameNumbers(X)
    if rem(size(X,2),3) == 0
        framenumbers = false;
    elseif rem(size(X,2),3) == 1
        % check to make sure it's not a zero column
        if X(2,1) - X(1,1) == 1
            framenumbers = true;
        else
            framenumbers = nan;
        end
    else
        framenumbers = nan;
    end

function [midpoints] = calcMidpoints(A, B, dist_thresh, max_pts)
% CalcNMidpoints calculate midpoints of 3D points
% Function generates a set of 3D points in between the two user-provided
% endpoints. The user can set the number of generated mid-point through a
% combination of a minimum distance threshold and/or a max number of
% points.
%
% Inputs
%   -A,B: endpoints, each of dimension 1x3 (xyz point)
%   -dist_thresh: minimum distance between points (can set to 0 if user
%   wants to use max_pts as the cutoff)
%   -max_pts: maximum number of points to calculated (can set to inf if
%   user wants to use a distance threshold as the cutoff)
%
% Output
%   -midpoints: calculated midpoints, where columns are x,y,z
%
%
%  written by J.D. Laurence-Chasen 4/13/2022
%

% make sure max-pts is an odd number
if max_pts < inf
    if rem(max_pts,2) == 0
        max_pts = max_pts + 1;
    end
end

d = norm(A-B) / 2; % initial distance

if d < dist_thresh
    warning("The (high) distance threshold for inter-locator bridges is preventing new locators from being created")
    midpoints = [];
else
    % calculate first midpoint
    
    pts = [A; mean([A;B]); B];
    npts = size(pts,1);
    d = norm(pts(1,:)-pts(2,:)) / 2;
    
    % add additional midpoints until threshold/number of points is reached
    while d > dist_thresh && (npts-2 + (npts - 3)) < max_pts
        
        newpts = [];
        
        for p = 1:npts-1
            
            p1 = pts(p,:);
            p2 = pts(p+1,:);
            p3 = mean([p1;p2]);
            
            newpts = [newpts(1:end-1,:); p1; p3; p2];
            
        end
        
        pts = newpts;
        
        % what would distance be with an additional division?
        d = norm(pts(1,:)-pts(2,:)) / 2;
        npts = size(pts,1);
    end
   midpoints = pts(2:end-1,:); 
end
