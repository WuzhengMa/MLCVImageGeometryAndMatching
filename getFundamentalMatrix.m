function [F, inliers] = getFundamentalMatrix(pts1, pts2, method, numTrials)
    integerClass = 'int32';
    outputClass = 'double';
    distanceThreshold = 0.01;
    conf = 0.99;

    threshold = cast(distanceThreshold, outputClass);
    nTrials = cast(numTrials, integerClass);
    numPts = cast(size(pts1, 1), integerClass);
    inliers = false(numPts, 1);
    homoPts1 = coder.nullcopy(zeros(3, numPts, outputClass));
    homoPts2 = coder.nullcopy(zeros(3, numPts, outputClass));
    homoPts1(3, :)   = 1;
    homoPts2(3, :)   = 1;

    homoPts1(1:2, :) = pts1';
    homoPts2(1:2, :) = pts2';

    if strcmp(method, 'Norm8Points')
        F = norm8Point(homoPts1, homoPts2, outputClass, integerClass);
        inliers = ones(1, numPts);
    elseif strcmp(method, 'RANSAC') 
        inliers(:) = ransac(outputClass, integerClass, homoPts1, homoPts2, numPts, nTrials,...
        threshold, conf);
        % Compute the fundamental matrix from the inliers
        F = norm8Point(homoPts1(:, inliers), homoPts2(:, inliers), outputClass,...
        integerClass);
    end
end

% Compute the distance of points according to a fundamental matrix.
function d = computeDistance(homoPts1, homoPts2, f)
pfp = (homoPts2' * f)';
pfp = pfp .* homoPts1;
d = sum(pfp, 1) .^ 2;

epl1 = f * homoPts1;
epl2 = f' * homoPts2;
d = d ./ (epl1(1,:).^2 + epl1(2,:).^2 + epl2(1,:).^2 + epl2(2,:).^2);

end

% Randomly select 8 points and compute the fundamental matrix.
function [d, f] = estTFormDistance(homoPts1, homoPts2, numPts,...
    outputClass, integerClass)

    indices = cast(randperm(numPts, 8), integerClass);
    f = norm8Point(homoPts1(:, indices), homoPts2(:, indices), outputClass,...
      integerClass);
    d = computeDistance(homoPts1, homoPts2, f);
end


% Function norm8Point computes the fundamental matrix using THE NORMALIZED
% 8-POINT ALGORITHM as described in page 281 of the following reference:
%   R. Hartley, A. Zisserman, "Multiple View Geometry in Computer Vision," 
%   Cambridge University Press, 2003. 
function f = norm8Point(homoPts1, homoPts2, outputClass, integerClass)
% Normalize the points
num = cast(size(homoPts1, 2), integerClass);
[homoPts1, t1] = vision.internal.normalizePoints(homoPts1, 2, outputClass);
[homoPts2, t2] = vision.internal.normalizePoints(homoPts2, 2, outputClass);

% Compute the constraint matrix
m = coder.nullcopy(zeros(num, 9, outputClass));
for idx = 1: num
  m(idx,:) = [...
    homoPts1(1,idx)*homoPts2(1,idx), homoPts1(2,idx)*homoPts2(1,idx), homoPts2(1,idx), ...
    homoPts1(1,idx)*homoPts2(2,idx), homoPts1(2,idx)*homoPts2(2,idx), homoPts2(2,idx), ...
                 homoPts1(1,idx),              homoPts1(2,idx), 1];
end

%Using svd
[~, ~, vm] = svd(m, 0);
f = reshape(vm(:, end), 3, 3)';

% Enforce rank-2 constraint
[u, s, v] = svd(f);
s(end) = 0;
f = u * s * v';

% Transform the fundamental matrix back to its original scale.
f = t2' * f * t1;

% Normalize the fundamental matrix.
f = f / norm(f);
if f(end) < 0
  f = -f;
end
end

 
function [inliers, nInliers] = findInliers(distance, numPts, threshold)
inliers = distance <= threshold;
nInliers = cast(sum(inliers), 'like', numPts);
end


 
% Update the number of trials based on the desired confidence and the
% inlier ratio.
function maxNTrials = updateNumTrials(oneOverNPts, logOneMinusConf, ...
  outputClass, integerClass, curNInliers, maxNTrials)

ratioOfInliers = cast(curNInliers, outputClass) * oneOverNPts;
if ratioOfInliers > cast(1, outputClass) - eps(outputClass)
  newNum = zeros(1, integerClass);
else
  ratio8 = ratioOfInliers^8;
  if ratio8 > eps(ones(1, outputClass))
    logOneMinusRatio8 = log(ones(1, outputClass) - ratio8);
    newNum = cast(ceil(logOneMinusConf / logOneMinusRatio8), integerClass);
  else
    newNum = intmax(integerClass);
  end
end

if maxNTrials > newNum
  maxNTrials = newNum;
end
end

function inliers = ransac(outputClass,...
    integerClass, homoPts1, homoPts2, numPts, nTrials, threshold, confidence)

    inliers = false(1, numPts);
    if numPts >= 8
      maxNTrials = nTrials;
      curNTrials = zeros(1, integerClass);
      bestNInliers = zeros(1, integerClass);
      logOneMinusConf = log(ones(1, outputClass) - confidence);
      oneOverNPts = ones(1, outputClass) / cast(numPts, outputClass);

      while curNTrials < maxNTrials
        d = estTFormDistance(homoPts1, homoPts2, numPts, outputClass,...
          integerClass);

        [curInliers, curNInliers] = findInliers(d, numPts, threshold);

        if bestNInliers < curNInliers
          bestNInliers = curNInliers;
          inliers = curInliers;

          % Update the number of trials
          maxNTrials = updateNumTrials(oneOverNPts, logOneMinusConf, ...
            outputClass, integerClass, curNInliers, maxNTrials);
        end
        curNTrials = curNTrials + 1;
      end
      
    else
      disp('Not enough points');
    end
end

