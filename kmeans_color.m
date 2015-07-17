%========================================================================%
%    Clusters regions using k-means                                      %
%                                                                        %
%    Author: Saket Navlakha                                              %
%    Date  : August 2010                                                 %
%                                                                        %
%========================================================================%


function L = kmeans_color(I,I_texture,L)
%KMEANS_COLOR   K-means clustering of regions.
%
% KMEANS_COLOR(I,L) computes a vector for each region corresponding to the
% region's normalized intensity histogram. Then applies k-means to cluster
% the vectors. Each cluster is colored using LABEL2RGB and the resulting
% color image is superimposed on top of the original image to produce a
% nicer-looking segmentation.
%
% >> I_texture = MR8fast(result.im)'; 
% Adding texture seems to help.

% ISBI:  kmeans_color(result.im,I_texture,post_process(Lx_Us.m1000,5));

%% Set constants.
%num_clusters = 13;   % # of k-means clusters.
num_clusters = 500;   % # of k-means clusters.
num_histogram_bins = 32;    % # of bins in the intensity histogram.
intensity_bins = 1:(256/num_histogram_bins):257; % actual bins.
alpha = 0.1;    % controls transparency. higher => more labels.
% was 0.2 for presentation one.


%% Compute region properties.
[L,num_regions] = bwlabel(L);
Regions = regionprops(L,I,'PixelValues','PixelIdxList','Area','MajorAxisLength','MinorAxisLength','Perimeter');
%Vectors = zeros(num_regions,num_histogram_bins+1);
Vectors = zeros(num_regions,num_histogram_bins+1+8); % texture.
%   Vectors = zeros(num_regions,4); % texture + shape

% TODO: see if texture can be used to differentiate between membranes and
% mitochondria.
 
for i=1:num_regions
    h = histc(Regions(i).PixelValues,intensity_bins); % edge centers.
    h = h./sum(h); % convert to probability distribution.
    
    % Texture.
    region_texture = I_texture(Regions(i).PixelIdxList,:);
    avg_region_texture = mean(region_texture,1);
    %size(h)
    if size(h,1) == 1
        h = h';
    end

    %Vectors(i,:) = h;
    Vectors(i,:) = [h' avg_region_texture];
    %Vectors(i,:) = [Regions(i).Area Regions(i).MajorAxisLength Regions(i).MinorAxisLength Regions(i).Perimeter];
    
    %if size(h,1)+size(avg_region_texture,2) ~= num_histogram_bins+9
    %    error('something')
    %end
end



num_regions


%% Apply k-means and re-number the label matrix.
C = kmeans(Vectors,num_clusters);

for i=1:num_regions
    L(Regions(i).PixelIdxList) = C(i);
end


L = remove_boundaries(L);


%% Superimpose the label matrix on top of the original image.
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure, imshow(I), hold on
himage = imshow(Lrgb);
set(himage, 'AlphaData', alpha);

end

% idea is to look at each boundary pixel in the given L. call the pixel to
% the left of it l, and right r. assuming L has been bwlabel'd correctly,
% check if C(l) == C(r). and if so, remove the boundary. do this
% recursively (maybe even call remove_boundaries).