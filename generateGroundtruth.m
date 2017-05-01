function [X1, Y1, X2, Y2] = generateGroundtruth(imageName1, imageName2, method, iter)

    patchSize = 32;
    HarrisThreshold = 500;
    HarrisRadius = 10;
    
    if size(size(imread(imageName1)),2) == 2
        imgExample1 = (imread(imageName1));
        imgExample2 = (imread(imageName2));
    else 
        imgExample1 = rgb2gray(imread(imageName1));
        imgExample2 = rgb2gray(imread(imageName2));
    end
    
    if(strcmp(method, 'Harris'))
        %Find interest points of images
        [x1, y1] = harrisDetector(imageName1, patchSize, HarrisThreshold, HarrisRadius);
        [x2, y2] = harrisDetector(imageName2, patchSize, HarrisThreshold, HarrisRadius);

        %Get color histogram descriptor
        colorHistogram = true; %Using the intensity at each pixel of the map is better than using color histogram
        descriptors1 = getDescriptors(imgExample1, x1, y1, patchSize, colorHistogram);
        descriptors2 = getDescriptors(imgExample2, x2, y2, patchSize, colorHistogram);

        %Match interest points by Harris
        [~, ~, x1, y1, x2, y2] = ...
            matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2, 'RANSAC');
    elseif(strcmp(method, 'SIFT'))
        [features1, descriptors1] = vl_sift(single(imgExample1), 'EdgeThresh', 3, 'PeakThresh', 10);
        [features2, descriptors2] = vl_sift(single(imgExample2), 'EdgeThresh', 3, 'PeakThresh', 10);
        
        [~, ~, x1, y1, x2, y2] = matchSIFTDescriptors(descriptors1, descriptors2,...
            features1(1,:), features1(2,:), features2(1,:), features2(2,:)); 
        
    end
    
    [h, inliers, HMstatus] = getHomographyMatrix(x1, y1, x2, y2, 'RANSAC', iter);
    
    X1 = x1(inliers);
    Y1 = y1(inliers);
    X2 = x2(inliers);
    Y2 = y2(inliers);
    
    subplot(1,2,1);
    imshow(imageName1);
    title('Interest points auto detection for image 1');
    hold on;
    %plot(features1(1,:), features1(2,:), 'rx');
    plot(Y1,X1,'rx');
    hold off; 

    subplot(1,2,2);
    imshow(imageName2);
    title('Interest auto detection for image 2');
    hold on;
    %plot(features2(1,:), features2(2,:), 'rx');
    plot(Y2,X2,'rx');
    hold off;
    
    text = sprintf('%d inliers out of %d interest points', size(X1,2), size(x1,2));
    disp(text);
end