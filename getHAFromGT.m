function [HA, HD] = getHAFromGT(H, GTx1, GTy1, GTx2, GTy2)
    %HD is computed as the sum of difference between transformed points and
    %the original points
    %HA is computed as the rate where the transformed points is identical
    %to the original points
    HA = 0;
    
    %Obtain the projection points from image 2 to image 1
    homoTransPoints = H * [GTx1;GTy1;ones(1,size(GTx1,2))];
    oneOverHomoZ=(1./homoTransPoints(3,:));
    oneOverHomoZ=[oneOverHomoZ; oneOverHomoZ; oneOverHomoZ];
    transPoints = oneOverHomoZ.*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates
    
    HD = abs([GTx2; GTy2] - transPoints(1:2,:));
    correctMapped = HD < ones(2, size(HD, 2));
    for i = 1:size(correctMapped,2)
        if correctMapped(1,i) && correctMapped(1,i)
            HA = HA + 1;
        end
    end
    HA = HA / size(correctMapped,2); %Homography Accuracy
    HD = sum(sum(HD))/size(GTx1,2);
    
end