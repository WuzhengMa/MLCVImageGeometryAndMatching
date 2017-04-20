function [e1, e2] = getEpipoles(F)
	[U,~,V] = svd(F,0);
	e1 = V(:,3);
    e1 = e1/e1(3,1);%Change back from homogenious coordinates
	e2 = U(:,3);
    e2 = e2/e2(3,1);%Change back from homogenious coordinates
end