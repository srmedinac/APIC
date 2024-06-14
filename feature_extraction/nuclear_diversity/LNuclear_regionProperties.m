%% borrow from Veta
function [P,N] = LNuclear_regionProperties(I, L)
if size(I,3)>1
    I=I(:,:,1);
end

I=double(I);

PADDING = 1.5;

props = {
    'Area'
    'Perimeter'
    'Eccentricity'
    'MajorAxisLength'
    'MinorAxisLength'
    'EquivDiameter'
    'Solidity'
    'Orientation'
    'Centroid'
    'WeightedCentroid'
    'FilledImage'
    'BoundingBox'
    };

% NOTE: the image is inverted
P = regionprops(L, 255-I, props);

% NOTE: the image is inverted
% extract separate images and label images (masks) for each region
[R M D] = labelmatrixToROI(255-I, L, PADDING, [], P);


N={};
% calculate additional properties
for k = 1:length(P)

    bw=bwboundaries(P(k).FilledImage);
    bw=bw{1,1};
    bbox=P(k).BoundingBox;
    N=[N;[bw(:,1)+bbox(2)-.5 bw(:,2)+bbox(1)-.5]];
    
    c = P(k).Centroid;
    w = P(k).WeightedCentroid;
    d = P(k).EquivDiameter;
    a = P(k).Area;
    l = P(k).Perimeter;
    
    region = R{k};
    regionMask = M{k};
    distance = D{k};
    
    pixelValues = region(regionMask);
    
    % circularity:
    circularity = 4 * pi * a / l^2;
    
    % elliptical deviation:
    % difference between the region and the elliptical approximation of the
    % region
    ellipse = distance <= 1;
    intersection = ellipse & regionMask;
    ellipticalDeviation = 1 - 2*sum(intersection(:)) / (P(k).Area + sum(ellipse(:)));
    
    % mass displacement:
    % distance betweend centroid and weighted centroid normalized with
    % equivalent diameter
    massDisplacement = sqrt(sum((c - w).^2))/d;
    
    % integrated intensity:
    integratedIntensity = sum(pixelValues);
    
    % mean intensity:
    meanIntensity = mean(pixelValues);
    
    % intensity deviation:
    intensityDeviation = std(pixelValues);
    
    % intensity range:
    intensityRange = prctile(pixelValues, 97.5)...
        - prctile(pixelValues, 2.5);
    
    % the inside boundary is defined as the residual between the image
    % region and its erosion with an isotropic strucutring element with
    % radius 1/8 of the equivalent diameter of the region
    se = strel('disk', round(d/8), 0);
    insideBoundary = xor(regionMask, imerode(regionMask, se));
    outsideBoundary = xor(regionMask, imdilate(regionMask, se));
    
    insideBoundaryValues = region(insideBoundary);
    outsideBoundaryValues = region(outsideBoundary);
    
    % inside boundary intensity statistics:
    meanInsideBoundaryIntensity = mean(insideBoundaryValues );
    insideBoundaryIntensityDeviation = std(insideBoundaryValues );
    insideBoundaryIntensityRange = prctile(insideBoundaryValues , 97.5) - prctile(insideBoundaryValues , 2.5);
    normalizedInsideBoundaryIntensity = meanInsideBoundaryIntensity / meanIntensity;
    
    % outside boundary intensuty statistics:
    meanOutsideBoundaryIntensity = mean(outsideBoundaryValues );
    outsideBoundaryIntensityDeviation = std(outsideBoundaryValues );
    outsideBoundaryIntensityRange = prctile(outsideBoundaryValues , 97.5) - prctile(outsideBoundaryValues , 2.5);
    normalizedOutsideBoundaryIntensity = meanOutsideBoundaryIntensity / meanIntensity;
    
    % boundary saliency:
    boundarySaliency = meanInsideBoundaryIntensity - meanOutsideBoundaryIntensity;
    normalizedBoundarySaliency = normalizedInsideBoundaryIntensity - normalizedOutsideBoundaryIntensity;
    
    % add to structure
    P(k).Circularity = circularity;
    P(k).EllipticalDeviation = ellipticalDeviation;
    P(k).MassDisplacement = massDisplacement;
    P(k).IntegratedIntensity =  integratedIntensity;
    P(k).MeanIntensity = meanIntensity;
    P(k).IntensityDeviation = intensityDeviation;
    P(k).IntensityRange = intensityRange;
    P(k).MeanInsideBoundaryIntensity = meanInsideBoundaryIntensity;
    P(k).InsideBoundaryIntensityDeviation = insideBoundaryIntensityDeviation;
    P(k).InsideBoundaryIntensityRange = insideBoundaryIntensityRange;
    P(k).NormalizedInsideBoundaryIntensity = normalizedInsideBoundaryIntensity;
    P(k).MeanOutsideBoundaryIntensity = meanOutsideBoundaryIntensity;
    P(k).OutsideBoundaryIntensityDeviation = outsideBoundaryIntensityDeviation;
    P(k).OutsideBoundaryIntensityRange = outsideBoundaryIntensityRange;
    P(k).NormalizedOutsideBoundaryIntensity = normalizedOutsideBoundaryIntensity;
    P(k).BoundarySaliency = boundarySaliency;
    P(k).NormalizedBoundarySaliency = normalizedBoundarySaliency;
    
end

end