function [HA, HD] = getHomoAccuracy(oriPoints, transPoints)
    %HD is computed as the sum of difference between transformed points and
    %the original points
    %HA is computed as the rate where the transformed points is identical
    %to the original points
    HA = 0;
    HD = sum(sum(abs(oriPoints - transPoints)));
    correctMapped = HD < [1;1];
    for i = 1:size(correctMapped,2)
        if correctMapped(1,i) && correctMapped(1,i)
            HA = HA + 1;
        end
    end
    HA = HA / size(correctMapped,2); %Homography Accuracy
end