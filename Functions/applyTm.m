function out = applyTm(point, mat)
% This function applies a 4x4 transformation matrix or a 3x3 rotation matrix
% to a point. Either the point or matrix can be static or moving. The
% function assumes rows = frames for the point, and z dim = frames for the
% matrix. Point format should always be xyz, and the code will add another
% parameter if needed (i.e. if its a 4x4 tranformation matrix). It will use mm2tm
% to convert RBT if coming from XMALab.
%
% Inputs - point = [x y z] where rows = frames
%        - mat = 4x4 transformation matrix or 3x3 rotation matrix *nframes
%
%
% Out = [x y z] * nframes
%
% Written by J.D. Laurence-Chasen, 2018


if size(mat,2) == 16
    mat = mm2tm(mat); % convert from maya matrix form
end

if size(mat,1) == 3
    tm = 0; % it's a rotation matrix
elseif size(mat,1) == 4
    tm = 1; % its a tranformation matrix
else 
    error('Check your matrix pls.')
end


% If point is static or it's just 1 frame
if size(point,1) == 1 && size(point,2) == 3
    
    point = point';
    nframes = size(mat,3);
    
    if tm == 1
        point = [point; 1];
        out = zeros(nframes,4);
    else
        out = zeros(nframes,3);
    end
    
    for i = 1:nframes
        out(i,:) = mat(:,:,i) * point;
    end
    
    
% If mat is static
elseif size(point,1) > 1 && size(point,2) == 3 && size(mat,3) == 1
    
    nframes = size(point,1);
    
    if tm == 1
        point(:,4) = 1;
        out = zeros(nframes,4);
    else 
        out = zeros(nframes,3);
    end
    
    for i = 1:nframes
        out(i,:) = mat * point(i,:)';
    end
    
% if both point and tm are dynamic
elseif size(point,1) > 1 && size(point,2) == 3 && size(mat,3) == size(point,1)
    
    nframes = size(point,1);
    
    if tm == 1
        point(:,4) = 1;
        out = zeros(nframes,4);
    else
        out = zeros(nframes,3);
    end
    
    for i = 1:nframes
        out(i,:) = mat(:,:,i) * point(i,:)';
    end
    
else
    error('Check input arguments for correct format')
end

% get rid of last column of 1s
if tm == 1
    out = out(:,1:3);
end

end

