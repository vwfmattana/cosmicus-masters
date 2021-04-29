clear;
%--------------------Paths-------------------------------------------------
path='F:\Ground Truth\Generate\Binarize\';
[filename, pathname] = ...
     uigetfile({'*.bmp;*.jpg;',...
 'Picture Files (*.bmp,*.jpg)';},'Select source image for processing');

%Input Files, and details thereof
%--------------------------------------------------------------------------
input=imread(strcat(strcat(pathname,'\'),filename) );
dimension=size(input);
Rdim=dimension(1);
Cdim=dimension(2);

fun3=@(block_struct) scalelines(block_struct.data), 'BorderSize', [5 5];
sl =blockproc( input, [Rdim,round(Cdim/18)], fun3);
imwrite(sl,strcat(path,'04-Scale Lines.bmp'),'bmp');