function Features = Generate_Features(I,L)
%GENERATE_FEATURES   Generates average intensity and texture features for
%each region.
%
% GENERATE_FEATURES(I,L) generates an NxK features matrix, where N is the
% number of regions in the label matrix L; and K is the 9 features for the
% region (avg intensity and avg texture for each of the 8 MR8 textons).
%

addPath('texture/');

num_regions = max(L(:));
Features = zeros(num_regions,9);

Regions = regionprops(L,I,'PixelValues','PixelIdxList'); % for intensity and texture features.
I_texture = MR8fast(I)';
    
for reg=1:num_regions
    
    % Compute average intensity for pixels in the region.
    avg_intensity = mean(Regions(reg).PixelValues);
    
    % Compute average texture along each dimension for pixels in the
    % region.
    region_texture = I_texture(Regions(reg).PixelIdxList,:); % # region pixels x 8.
    avg_texture = mean(region_texture,1);
    
    Features(reg,:) = [avg_intensity avg_texture];
end