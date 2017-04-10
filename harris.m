function [cim, r, c] = harris(im, thresh, radius)
    
    error(nargchk(2,5,nargin));
    
    dx = [-1 0 1; -1 0 1; -1 0 1]; % Derivative masks
    dy = dx';
    
    Ix = conv2(im, dx, 'same');    % Image derivatives
    Iy = conv2(im, dy, 'same');    

    % Generate Gaussian filter of size 6*sigma (+/- 3sigma) and of
    % minimum size 1x1.
    %g = fspecial('gaussian',max(1,fix(6*sigma)), sigma);
    temp = [0.03 0.105 0.222 0.286 0.222 0.105 0.03];
    g = temp'*temp;
    
    Ix2 = conv2(Ix.^2, g, 'same'); % Smoothed squared image derivatives
    Iy2 = conv2(Iy.^2, g, 'same');
    Ixy = conv2(Ix.*Iy, g, 'same');
    
    cim = ((Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps)); % Harris corner measure

    % Alternate Harris corner measure used by some.  Suggested that
    %k=0.04; %- I find this a bit arbitrary and unsatisfactory.
    %cim = (Ix2.*Iy2 - Ixy.^2) - k*(Ix2 + Iy2).^2; 

	
	% Extract local maxima by performing a grey scale morphological
	% dilation and then finding points in the corner strength image that
	% match the dilated image and are also greater than the threshold.
	sze = 2*radius+1;                   % Size of mask.
	mx = ordfilt2(cim,sze^2,ones(sze)); % Grey-scale dilate.
	cimmx = (cim==mx)&(cim>thresh);       % Find maxima.
	
	[r,c] = find(cimmx);                  % Find row,col coords.
    end