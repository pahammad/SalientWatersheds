% code to plot edges in red on top of original image
% Parvez Ahammad, HHMI/JFRC

function I_s = plotGreenEdges(originalIm, edgeIm, plot_flag)

%originalIm = adjustImageContrast(originalIm);

I_s = uint8(repmat(originalIm, [1 1 3]));
edgeIm = edgeIm>0;
G = I_s(:,:,2);
G(edgeIm) = 255;
I_s(:,:,2) = G;
if plot_flag == 1
    figure, imshow(I_s)
end
