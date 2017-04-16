function [sSizeDesp, lSizeDesp, SX, SY, LX, LY] = keepEqualNumPoints(sSizeDescriptors, lSizeDescriptors, sx, sy, lx, ly, inlinerFactor)
    %Use NN search to find equal pairs of closest interest points
    sSizeDesp = [];
    lSizeDesp = [];
    SX = [];
    SY = [];
    LX = [];
    LY = [];
    [nearestIndex, distance] = knnsearch(lSizeDescriptors, sSizeDescriptors, 'Distance', 'cityblock');
    inlinerFlag = distance < mean(distance)*inlinerFactor;
    for i = 1:size(inlinerFlag)
        if inlinerFlag(i)
            sSizeDesp = sSizeDescriptors;
            lSizeDesp = [lSizeDesp; lSizeDescriptors(nearestIndex(i),:)];
            SX = [SX, sx(i,:)'];
            SY = [SY, sy(i,:)'];
            LX = [LX, lx(nearestIndex(i), :)];
            LY = [LY, ly(nearestIndex(i), :)];
        end
    end
end