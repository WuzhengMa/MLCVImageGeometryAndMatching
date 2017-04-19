init;
patchSize = 32;

%% Q2.1 Find interest points on origin/resized image
imageName1 = 'HGFigures/DSC02715.ppm';
imageName1ResizedHalfColor = 'HGFigures/DSC02715ResizedHalf.ppm';
%imageName1ResizedHalfColor = (imread(imageName1));
%imageName1ResizedHalfColor = imresize(imageName1ResizedHalfColor, 0.5);
%imwrite(imageName1ResizedHalfColor,'HGFigures/DSC02715ResizedHalf.ppm','ppm')

if size(size(imread(imageName1)),2) == 2
    imgExample1 = (imread(imageName1));
    imgExample1ResizedHalf = (imread(imageName1ResizedHalfColor));
else 
    imgExample1 = rgb2gray(imread(imageName1));
    imgExample1ResizedHalf = rgb2gray(imread(imageName1ResizedHalfColor));
end

%imshow(imgExample1);

%Find interest points of images
[x1, y1] = harrisDetector(imageName1, patchSize);
[x1ResizedHalf, y1ResizedHalf] = harrisDetector(imageName1ResizedHalfColor, patchSize);

%Get color histogram descriptor
colorHistogram = true; %Using the intensity at each pixel of the map is better than using color histogram
descriptors1 = getDescriptors(imgExample1, x1, y1, patchSize, colorHistogram);
descriptors1ResizedHalf = getDescriptors(imgExample1ResizedHalf, x1ResizedHalf, y1ResizedHalf, patchSize, colorHistogram);


%Show interest points
subplot(2,2,1);
imshow(imageName1);
title('Interest points detection for original image');
hold on;
%plot(features1(1,:), features1(2,:), 'rx');
plot(y1,x1,'rx');
hold off; 

subplot(2,2,2);
imshow(imageName1ResizedHalfColor);
title('Interest points detection for resized image');
hold on;
%plot(features2(1,:), features2(2,:), 'rx');
plot(y1ResizedHalf,x1ResizedHalf,'rx');
hold off;


%Match interest points
[descriptors1, descriptors1ResizedHalf, x1, y1, x1ResizedHalf, y1ResizedHalf] = ...
    matchDescriptorSize(descriptors1, descriptors1ResizedHalf, x1, y1, x1ResizedHalf, y1ResizedHalf, 'Norm8Points');

%Show interest points
subplot(2,2,3);
imshow(imageName1);
title('Interest points detection for original image');
hold on;
%plot(features1(1,:), features1(2,:), 'rx');
plot(y1,x1,'rx');
hold off; 

subplot(2,2,4);
imshow(imageName1ResizedHalfColor);
title('Interest points detection for resized image');
hold on;
%plot(features2(1,:), features2(2,:), 'rx');
plot(y1ResizedHalf,x1ResizedHalf,'rx');
hold off;

%Calculate homography
h = getHomographyMatrix(x1, y1, x1ResizedHalf, y1ResizedHalf, 'RANSAC', 2000);

%Obtain the projection points from image 2 to image 1
homoTransPoints = h\[x1ResizedHalf;y1ResizedHalf;ones(1,size(x1ResizedHalf,2))];
oneOverHomoZ=(1./homoTransPoints(3,:));
oneOverHomoZ=[oneOverHomoZ; oneOverHomoZ; oneOverHomoZ];
transPoints = oneOverHomoZ.*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates

%Calculate the homography accuracy HA
[HA, HD] = getHomoAccuracy([x1; y1], transPoints([1,2], :))
