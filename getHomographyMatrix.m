function [res, inliers, status] = getHomographyMatrix(x1, y1, x2, y2, method, nTrials)
    %[x2, y2] = rearrangePoints(nearestIndex, x2, y2);
    confidence = 99;
    threshold = 0.01;
    
    statusCode = struct('NoError',           int32(0),...
                    'NotEnoughPts',      int32(1),...
                    'NotEnoughInliers',  int32(2));
   
    n = size(x1,2);
    if(n >= 4)
        if strcmp(method, 'default')

            rows0 = zeros(3, n);
            rowsXY = -[x1; y1; ones(1,n)];
            hx = [rowsXY; rows0; x2.*x1; x2.*y1; x2];
            hy = [rows0; rowsXY; y2.*x1; y2.*y1; y2];
            h = [hx hy];
            [U, ~, V] = svd(h');
            res = (reshape(V(:,9), 3, 3)).';
            res = res./res(3,3);
            inliers(:) = true(1, n);
            status = statusCode.NoError;
        else
            if(size(x1,2) > size(x2,2))
                nPts = size(x2,2);
            else
                nPts = size(x1,2);
            end
            inliers = false(nPts, 1);
            maxNTrials = nTrials;
            curNTrials = 0;
            bestNInliers = 0;
            logOneMinusConf = log(1 - confidence);
            oneOverNPts = 1 / nPts;

              while curNTrials < maxNTrials
                [d, H] = estHomoMatrix(x1, y1, x2, y2);

                [curInliers, curNInliers] = findInliers(d, threshold);

                if bestNInliers < curNInliers
                    bestNInliers = curNInliers;
                    res = H;
                    inliers = curInliers;
                    % Update the number of trials
                    %maxNTrials = updateNumTrials(oneOverNPts, logOneMinusConf, curNInliers, maxNTrials);
                end
                curNTrials = curNTrials + 1;
              end

              if bestNInliers >= 4
                status = statusCode.NoError;
              else
                status = statusCode.NotEnoughInliers;
              end
              res = res./res(3,3);
        end
    else
        status = statusCode.NotEnoughInliers;
        inliers = false(n, 1);
        res=zeros(3,3);
    end


%========================================================================== 
% Randomly select 4 points and compute the homography matrix.
%========================================================================== 
function [d, H] = estHomoMatrix(x1, y1, x2, y2)

pts = datasample([x1; y1; x2; y2;] ,4 , 2,  'Replace', false);
H = norm4Point(pts(1,:), pts(2,:), pts(3,:), pts(4,:));
d = calHADistance(x1, y1, x2, y2, H);

function res = norm4Point(x1, y1, x2, y2)
    
rows0 = zeros(3, 4);
rowsXY = -[x1; y1; ones(1,4)];
hx = [rowsXY; rows0; x2.*x1; x2.*y1; x2];
hy = [rows0; rowsXY; y2.*x1; y2.*y1; y2];
h = [hx hy];
[U, ~, V] = svd(h');
res = (reshape(V(:,9), 3, 3)).';

function d = calHADistance(x1, y1, x2, y2, H)

%Obtain the projection points from image 2 to image 1
homoTransPoints = H\[x2;y2;ones(1,size(x2,2))];
oneOverHomoZ=(1./homoTransPoints(3,:));
oneOverHomoZ=[oneOverHomoZ; oneOverHomoZ; oneOverHomoZ];
transPoints = oneOverHomoZ.*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates

%Calculate the homography accuracy HA
disM = [x1; y1] - transPoints([1,2], :);
d = hypot(disM(1,:), disM(2,:));

function [curInliers, curNInliers] = findInliers(d, threshold)

idx = 1;
curNInliers = 0;
curInliers = false(size(d,2), 1);
while(idx <= size(d,2))
    if(d(idx) <= threshold)
        curNInliers = curNInliers + 1;
        curInliers(idx) = true;
    end
    idx = idx + 1;
end

%========================================================================== 
% Update the number of trials based on the desired confidence and the
% inlier ratio.
% Referenced from Matlab estimateFundamentalMatrix.m function
%========================================================================== 
function maxNTrials = updateNumTrials(oneOverNPts, logOneMinusConf, curNInliers, maxNTrials)

ratioOfInliers = curNInliers * oneOverNPts;
if ratioOfInliers > 1 - eps
  newNum = 0;
else
  ratio8 = ratioOfInliers^8;
  if ratio8 > eps
    logOneMinusRatio8 = log(1 - ratio8);
    newNum = ceil(logOneMinusConf / logOneMinusRatio8);
  else
    newNum = intmax;
  end
end

if maxNTrials > newNum
  maxNTrials = newNum;
end
