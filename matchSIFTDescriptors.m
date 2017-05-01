function [Desc1, Desc2, X1, Y1, X2, Y2] = matchSIFTDescriptors(descriptors1, descriptors2, x1, y1, x2, y2) 
    if(size(x1) <= size(x2))
        matches = vl_ubcmatch(descriptors1, descriptors2);
        Desc1 = descriptors1(matches(1,:));
        Desc2 = descriptors2(matches(2,:));
        X1 = x1(matches(1,:));
        Y1 = y1(matches(1,:));
        X2 = x2(matches(2,:));
        Y2 = y2(matches(2,:));
    else
        matches = vl_ubcmatch(descriptors2, descriptors1);
        Desc2 = descriptors2(matches(1,:));
        Desc1 = descriptors1(matches(2,:));
        X2 = x2(matches(1,:));
        Y2 = y2(matches(1,:));
        X1 = x1(matches(2,:));
        Y1 = y1(matches(2,:));
    end
end