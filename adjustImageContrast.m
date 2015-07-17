% code to plot edges in red on top of original image
% Parvez Ahammad, HHMI/JFRC

function adjIm = adjustImageContrast(Im)

nChannels = size(Im,3);

for i=1:nChannels
    foo = Im(:,:,i);
    curMin = min(foo(:));
    curMax = max(foo(:));
    adjIm(:,:,i) = round((foo-curMin)*255/(curMax-curMin));
    clear foo;
end
