function [x, y] = getManualInterestPoints(imageName)
    %Manually get coordinates of interest points
    imgExample = imread(imageName);
    imshow(imgExample);
    [y, x] = ginput;
    x = fix(x);
    y = fix(y);

end