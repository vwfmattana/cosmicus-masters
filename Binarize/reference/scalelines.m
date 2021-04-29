% =========================================================================
% Vincent Mattana 21128707 2016 MSc IT/Computer Science
% Scalelines
% =========================================================================
clear;
[filename, pathname] = ...
     uigetfile({'*.bmp;*.jpg;',...
 'Picture Files (*.bmp,*.jpg)';},'Select source image for processing');
input=imread(strcat(strcat(pathname,'\'),filename) );

%--------------------Paths-------------------------------------------------
path='F:\LosTheRed\Computer Science\Masters\cosmicus\2016\Artefact\Ground Truth\Generate\Scalelines\'; 


d=3;
dimension=size(input);
ndim=dimension(1);
mdim=dimension(2);

g_counter=zeros([ndim mdim]);
g_counter=padarray(g_counter, [d d], 0);
g_counter=double(g_counter);

%--------------------Approximation-------------------------------------
for m=1+d:mdim-1-d,
prev=input(1+d,m);
counter=0;
    for n=1+d:ndim-1-d,
        counter=counter+1;
        record(n,m)=counter;
        if (input(n,m)~=prev)
            if (prev==1)
                g_counter(n-1,m)=counter;
            else 
                g_counter(n,m)=counter;
            end
            prev=input(n,m);    
            counter=0;
        end
    end
end
% figure, imshow(g_counter/64);figure(gcf);
% image shows scaleline positions + weight,
% where weight=the distance from upper(previous) scale line
% *weight=1,2,3...  ;probably the bottom of some scale line
%       =    +-12   ;probably the top of a scale line
%       =   16+     ;probably scaleline,below some missing scaleline(s)

for m=1+d:mdim+2*d,
     scaleline_count_per_column(m-d)=nnz(g_counter(:,m));
end     %nnz=number-of-nonzero elements,
% count the number of possible scalelines represented by each column
% low count = column at position of hourmarker ?
% medium count = ...at pos of garbage/unreadable piece
% high count = should contain information on every scale line present

candidate_scaleline_colnr = find(scaleline_count_per_column==max(scaleline_count_per_column),1, 'first');
% find("where some array satisfies logic"  ,  "how many items?"  , "which element to return if multiples found" )
% returns the column number where most scalelines dwell
scaleline_imageaprox(:,1)=g_counter(:,candidate_scaleline_colnr);
% numel(scaleline_imageaprox)
% nnz(scaleline_imageaprox)
scaleline_imageaprox = imresize(scaleline_imageaprox, [ndim mdim], 'nearest');

scaleline_imageaprox=im2bw(scaleline_imageaprox, 0.1);
% create full size image of scalelines by using the candidate as template
% binarize scale line image, since counter left elements > [0:1]
imwrite(scaleline_imageaprox,strcat(path,'scaleline_imageaprox.bmp'),'bmp');

sl_fillaprox=imdilate(scaleline_imageaprox,strel('line',4,90));
sl_fillaprox=imdilate(sl_fillaprox,strel('line',4,90));
sl_fillaprox=imerode(sl_fillaprox,strel('line',4,90));
% fill gaps between scaleline markers(from g_counter)



sl_skelaprox=imerode(sl_fillaprox,strel('line',2,90));
sl_skelaprox=imerode(sl_skelaprox,strel('line',2,90));
sl_skelaprox=bwmorph(sl_skelaprox,'thin',5);
imwrite(sl_skelaprox,strcat(path,'sl_skelaprox.bmp'),'bmp');


num_of_SL=nnz(sl_skelaprox(:,candidate_scaleline_colnr));
c_last_SL=0;
p_between_sl=round(ndim/num_of_SL);

for n=1:ndim,       %-----add missing SL based on average SL seperation
    c_last_SL=c_last_SL+1;
    if (sl_skelaprox(n,candidate_scaleline_colnr)>0)
        c_last_SL=0;
    end
    
    if (c_last_SL>(p_between_sl)*1.5)
        sl_skelaprox(n-round(c_last_SL*0.333),:)=1;
        c_last_SL=0;
    end
end
imwrite(sl_skelaprox,strcat(path,'90deg 05_skelaprox.bmp'),'bmp');