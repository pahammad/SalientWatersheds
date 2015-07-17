% code to plot edges in red on top of original image
% Parvez Ahammad, HHMI/JFRC

function I_s = plotRedGreenEdges(originalIm, edgeIm1, edgeIm2, plot_flag)

%originalIm = adjustImageContrast(originalIm);

I_s = uint8(repmat(originalIm, [1 1 3]));
edgeIm1 = edgeIm1>0;
edgeIm2 = edgeIm2>0;
R = I_s(:,:,1);
R(edgeIm1) = 255;
I_s(:,:,1) = R;

G = I_s(:,:,2);
G(edgeIm2) = 255;
I_s(:,:,2) = G;

if plot_flag == 1
    figure, imshow(I_s)
end
