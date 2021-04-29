function output=datalines(input)
%°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º° 
%Function to find scale lines on historic cosmic ray recordings, using
%average intensity and Hough transform.
%Format:
%               output=datalines(input)
%
%Input:
%               input: Unprocessed historic cosmic ray recording.
%
%Output:
%               output: An image containing only the data lines.
%
%Usage example:
%               data = datalines(binarized);

%Data lines
%--------------------------------------------------------------------------
path='F:\Ground Truth\Generate\Binarize\';
dimension=size(input);
Rdim=dimension(1);
Cdim=dimension(2);

    data=input(:,:,3);
    for y=1:Cdim,
        for x=1:Rdim,
            if ((input(x,y,1)==1)||(input(x,y,2)==1))
                data(x,y)=0;
            end
        end
    end
    
    data= bwareaopen(data,25);
    data=edge(data,'canny');
    data= bwareaopen(data,30);
    data=imdilate(data,strel('octagon', 6));
    data=bwmorph(data,'skel',Inf);
    data=bwmorph(data, 'spur', 200);
    data= bwareaopen(data,30);
    data=imdilate(data,strel('octagon', 6));
    CC = bwconncomp(data);
    
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    data=zeros(size(data));
    data(CC.PixelIdxList{idx}) = 1;
    
    output=data;
end