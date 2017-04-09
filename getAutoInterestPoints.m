%Harris interest point detector
function [cim, x, y] = getAutoInterestPoints(imgExample, thresh, radius) 

[cim, x, y] = harris(imgExample, thresh, radius);

end