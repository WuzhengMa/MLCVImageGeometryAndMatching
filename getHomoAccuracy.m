function [HA, HD] = getHomoAccuracy(oriPoints, transPoints)
    %HD is computed as the sum of difference between transformed points and
    %the original points
    %HA is computed as the rate where the transformed points is identical
    %to the original points
    HA = 0;
    HD = abs(oriPoints - transPoints);
    correctMapped = HD < ones(2, size(HD, 2));
    for i = 1:size(correctMapped,2)
        if correctMapped(1,i) && correctMapped(1,i)
            HA = HA + 1;
        end
    end
    HA = HA / size(correctMapped,2); %Homography Accuracy
    HD = sum(sum(HD));
    
end