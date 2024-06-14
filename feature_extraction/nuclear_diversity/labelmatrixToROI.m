
function [R M D] = labelmatrixToROI(I, L, padding, interpType, P)
% labelmatrixToROI: Extracts regions from an image given with a labelmatrix
% as separate images such that the major and minor axes of the regions are
% aligned with the image axes.

if  ~exist('padding', 'var') || isempty(padding)
    padding = 1;
end

if  ~exist('interpType', 'var') || isempty(interpType)
    interpType = 'linear';
end

% the properties matrix can be given as an input if already computed
if  ~exist('P', 'var') || isempty(P)
    props = {
        'Centroid'
        'Orientation'
        'MajorAxisLength'
        'MinorAxisLength'
        };
    
    P = regionprops(L, props);
end

R = cell(length(P), 1);
M = cell(length(P), 1);
D = cell(length(P), 1);

iI = griddedInterpolant(I, interpType);
iL = griddedInterpolant(L, 'nearest');

warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

for k = 1:length(P)
    
    % location, size and orientation of the image
    c = P(k).Centroid;
    a = P(k).MajorAxisLength;
    b = P(k).MinorAxisLength;
    t = P(k).Orientation;
    
    w = round(a * padding);
    h = round(b * padding);
    phi = degToRad(t);
    
    % form a sampling grid
    x = -(h/2-0.5):(h/2-0.5);
    y = -(w/2-0.5):(w/2-0.5);
    
    [X Y] = ndgrid(x, y);
    
    % rotate and center
    A = rotationTranslationMatrix(phi, c);
    T = maketform('affine', A);
    
    [V U] = tformfwd(T, Y, X);
    
    % sample image and label matrix
    %region = interp2(I, U, V, interpType, 0);
    %regionMask = interp2(L, U, V, '*nearest', 0) == k;
    region = iI(U, V);
    regionMask = iL(U, V)==k;
    
    % NOTE: The sampling of the region mask assumes propper region
    % ordering: k-th region has label k.
    
    % normalized distance from center
    distance = sqrt((2*X/a).^2 + (2*Y/b).^2);
    
    region(isnan(region)) = 0;
    regionMask(isnan(regionMask)) = false;
    
    R{k} = region;
    M{k} = regionMask;
    D{k} = distance;
    
end

warning('on', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

end