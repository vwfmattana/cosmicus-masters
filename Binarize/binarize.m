function output=binarize(input)
%°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º° 
%Function to binarize historic cosmic ray recordings, using
%
%Format:
%               output=binarize(input)
%
%Input:
%               input: Unprocessed historic cosmic ray recording.
%
%Output:
%               output: A binarized image.
%
%Usage example:
%               binarized_0 = binarize(input);

%Scale lines
%--------------------------------------------------------------------------
path='F:\Ground Truth\Generate\Binarize\';
dimension=size(input);
Rdim=dimension(1);
Cdim=dimension(2);

thresh=graythresh(input);
binarized=im2bw(input,thresh);

for k=1:200
    whitecount = sum(sum(binarized == 1))/(Rdim*Cdim);
    if (whitecount > 0.20) && (thresh<1) && (thresh>0) 
        thresh = thresh+0.004;
    end
    if (whitecount < 0.15) && (thresh<1) && (thresh>0) 
        thresh = thresh-0.002;
    end
    if (thresh<1) && (thresh>0) 
        binarized = im2bw(input, thresh);
    end
    
end
    %---------------Removing Noise---------------------------
%     binarized=imclose(binarized,strel('rectangle',[2 2]));
    binarized = bwareaopen(binarized, 25);
%     binarized=imclose(binarized,strel('rectangle',[2 2]));

output=binarized;
end