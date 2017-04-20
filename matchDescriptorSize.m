function [desp1, desp2, X1, Y1, X2, Y2] = matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2, method)
    %When the method is 'Norm8Points', only 8 pairs of interest points will 
    %be used, hence use NN method to find 8 pairs of closest interest points
    %When the method is 'RANSAC', equal pairs of interest points will be
    %used, hence use NN method to find equal pairs of closest interest
    %points
    inlinerFactor = 1.8;
    if size(descriptors2, 1) >= size(descriptors1, 1)
        if strcmp(method, 'Norm8Points')
            [desp1, desp2, X1, Y1, X2, Y2] = keepEightNearestPoints(descriptors1, descriptors2, x1, y1, x2, y2, inlinerFactor);
        elseif strcmp(method, 'RANSAC')
            [desp1, desp2, X1, Y1, X2, Y2] = keepEqualNumPoints(descriptors1, descriptors2, x1, y1, x2, y2, inlinerFactor);
        elseif strcmp(method, 'SortedDist')
            [desp1, desp2, X1, Y1, X2, Y2] = keepEqualNumPointsSortedDist(descriptors1, descriptors2, x1, y1, x2, y2, inlinerFactor);
        end
    elseif size(descriptors2, 1) < size(descriptors1, 1)
        if strcmp(method, 'Norm8Points')
            [desp2, desp1, X2, Y2, X1, Y1] = keepEightNearestPoints(descriptors2, descriptors1, x2, y2, x1, y1, inlinerFactor);
        elseif strcmp(method, 'RANSAC')
            [desp2, desp1, X2, Y2, X1, Y1] = keepEqualNumPoints(descriptors2, descriptors1, x2, y2, x1, y1, inlinerFactor);
        elseif strcmp(method, 'SortedDist')
            [desp2, desp1, X2, Y2, X1, Y1] = keepEqualNumPointsSortedDist(descriptors2, descriptors1, x2, y2, x1, y1, inlinerFactor);
        end
    end
end