function [tblout] = generateLocatorFile(locatorsCTspace,rbt,fileout,frames)
%generateLocatorFile Generate locator XYZ file for dynamicEndocastByBone
%
% % NOTE: THIS FUNCTION IS STILL IN BETA, USE AT OWN RISK 
%
% locatorsCTspace: CTexp file path (.csv) OR cell array where each cell has 
% rbt: rbt matrix, from import_rbt() xromm-tools function
% fileout: character string of output filename
% frames: (optional) vector of frame numbers for final file leave empty []
% for default
%
% Importantly, order of RBTs must be the same--i.e., locator positions in
% cell 1 must correspond to first rbt. Function will try to automatically
% detect diff locators based on first 5 characters in name
%
%
% % Written by J.D. Laurence-Chasen 3/8/2022

try
    tbl = readtable(locatorsCTspace);
    isptsfile = true;
catch
    isptsfile = false;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Import CTexp file and read column headers % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isptsfile % if its a file, separate by column header prefix
    
    % first check to make sure column headers are characters
    try
        csvread(locatorsCTspace,0,0); % this will fail if there are headers
        error('No column headers in locator file')
    catch
    end
    
    varnames = tbl.Properties.VariableNames; % get column names
    XYZfile = csvread(locatorsCTspace,1,0);
    % let's also make sure it's CT exp format (i.e. 1 row of data)
    if size(XYZfile,1) ~= 1
         error('Hmmm, locator CT space file should only have 1 row of data')
    end
 
    croppedvars = varnames';
    for v = 1:length(varnames)
        croppedvars{v,1} = varnames{v}(1:end-2); % crop off _x/_y/_z (diff than maya version _tx _ty _tz)
    end
    varnamesnew = croppedvars;
    %varnamesnew = unique(croppedvars,'stable'); % get rid of duplicates
    % find numbers at end of cell
    nums = regexp(varnamesnew,'\d*','Match');
    varnames_nonums = cell(length(nums),1); % This is the new list
    for v = 1:length(nums)
        if ~isempty(nums{v})
            tmpv = varnamesnew{v}(1:end-length(nums{v}{1}));
            varnames_nonums{v,1} = tmpv;
        end
    end

    bones = unique(varnames_nonums,'stable');
    
    % now we have bones, let's extract xyz positions and place into
    % individual cells    
    xyz_by_bone = cell(length(bones),1);
    for b = 1:length(bones)
        xyz_by_bone{b} = XYZfile(1,strcmp(varnames_nonums,bones{b}));
    end
else
    error('For now, input needs to be a CTexp CSV file of locators')
end

% Now let's import/check RBT formatting
try
    if size(rbt,4) == length(bones)
        disp(['Applying RBTs to locators on ' num2str(length(bones)) ' bones'])
    else
        error('Hmmm...dectected a different number of bones in locator file and rbt array')
    end
catch
    error('rbt needs to be a 4-D array generated with import_rbt()')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  % % Animate locators % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
animatedlocs = [];
for b = 1:length(bones) % for each bone
    
    nlocs = size(xyz_by_bone{b},2) / 3;
    
    for loc = 1:nlocs % for each locator
        
        cols = (loc-1)*3+1:(loc-1)*3+3;
        newpts = applyTm(xyz_by_bone{b}(cols),rbt(:,:,:,b)); % animate
        animatedlocs = [animatedlocs newpts]; % add to main array
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % Construct output table % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(frames)
frames = 1:size(rbt,3);
frames = frames';
else
    if size(frames,1) == 1
        frames = frames';
    end
end

tblout = array2table([frames animatedlocs], 'VariableNames',[{'frame'} varnames]);
writetable(tblout,fileout);


end


