function [x,y] = missile_centroid(image)
%
% Centroid.m computes the centroid of the image
%
%		inputs:
%			image     =  image data.
%
%
    [N,M] = size(image);
    dx = 1;
    y = dx*(1:N)';
    yy = y*ones(1,M);
    
    x	 = 1:M;
    xx = ones(N,1)*x;
    
    Wt = sum(sum(image));
    
    if Wt == 0
        x = M/2;
        y = N/2;
    else
        x = sum(sum(xx.*image))/Wt;
        y = sum(sum(yy.*image))/Wt;
    end

   
   return     