('***************START*******************')
path='F:\Ground Truth\Generate\Binarize\';
[filename, pathname] = ...
    uigetfile({'*.bmp;*.jpg;',...
    'Picture Files (*.bmp,*.jpg)';},'Select ORIGINAL image for processing');
original=imread(strcat(strcat(pathname,'\'),filename) );

path='F:\Ground Truth\Generate\Binarize\';
[filename, pathname] = ...
    uigetfile({'*.bmp;*.jpg;',...
    'Picture Files (*.bmp,*.jpg)';},'Select USER image for processing');
test=imread(strcat(strcat(pathname,'\'),filename) );

MS_Err = immse(original, test)
Cdim=size(original,1);
Rdim=size(original,2);

correct_assignment=0;
false_negative =0;
false_positive =0;
true_positive=0;
true_negative=0;
for ii= 1:Cdim,
    for jj= 1:Rdim,
        if (test(ii, jj) == original(ii, jj))
            correct_assignment=correct_assignment+1;
            if (test(ii, jj) ~= 0 && original(ii, jj) ~= 0)
                true_positive=true_positive+1;
            else
                true_negative=true_negative+1;
            end
        else
            if (test(ii, jj) == 0 && original(ii, jj) ~= 0)
                false_negative=false_negative+1;
            else
                false_positive=false_positive+1;
            end
        end
    end
end
total =Rdim*Cdim
correct_assignment=correct_assignment/total
disp('***************************************')
true_positive
true_negative
false_negative
false_positive
disp('***************************************')
precision=true_positive/(true_positive+false_positive)
recall=true_positive/(true_positive+false_negative)
accuracy=(true_positive+true_negative)/(true_positive+true_negative+false_positive+false_negative)
F_measure=2*((precision*recall)/(precision+recall))
