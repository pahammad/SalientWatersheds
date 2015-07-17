%========================================================================%
%    VALIDATION MEASURES                                                 %
%                                                                        %
%    Author: Saket Navlakha                                              %
%    Date  : August 2010                                                 %
%                                                                        %
%========================================================================%


function [L_gt,L_mistakes] = validation(GT_image,L_alg,connsize)
%VALIDATION calls the validation functions.
% See boundary_dist, asym_dist, partition_dist, and other_indices.
%

tic;

%% set flags.
boundary_dist_flag = 1;
asym_dist_flag = 1;
partition_dist_flag = 1;
other_indices_flag = 0;

red_pixel_indicator = 100; % below this and it's white.
L_mistakes = 0;


%% Convert ground truth image to label matrix.
fprintf('\tConverting GT colored image to label matrix...\n')

[num_rows, num_cols, ~] = size(GT_image);
L_gt = ones(num_rows,num_cols);

for i=1:num_rows
    for j=1:num_cols

        % check if it's a red pixel; if so, mark it a 0.
         if (GT_image(i,j,2) < red_pixel_indicator && GT_image(i,j,3) < red_pixel_indicator)
             L_gt(i,j) = 0;
         end
%          if (I(i,j,2) < 200 && I(i,j,3) < 200 )
%              L1(i,j) = 0; % boundary.
%          end
    end    
end

L_gt = bwlabel(L_gt,4);

% thin the boundaries -- also look at non-maxima suppression.
% TODO: are more regions added after this operation? YES.
L_gt = imcomplement(L_gt>0);
L_gt = imcomplement(bwmorph(L_gt,'skel',Inf));


%% call functions.

if boundary_dist_flag
    fprintf('\tComputing boundary distance...\n')
    [L_mistakes,correct,incorrect] = boundary_dist(L_gt,L_alg);

    fprintf('\t\tAccuracy = %.4f\n',correct/(correct+incorrect))
end


if (asym_dist_flag || partition_dist_flag)
    fprintf('\tComputing M...\n')
    M = compute_m(L_gt,L_alg,connsize);

    if asym_dist_flag
        fprintf('\tComputing asymmetric distance...\n')
        [n2,ad] = asym_dist(M,L_gt,L_alg);
        %fprintf('\t\tPercent of pixels to delete from L2 = %.4f (of %i)\n',ad/n2,n2)
        fprintf('\t\tPercent of pixels that match = %.4f (of %i)\n',1-ad/n2,n2)
    end

    if partition_dist_flag
        fprintf('\tComputing partition distance...\n')
        [n1,n2,pd] = partition_dist(M,L_gt,L_alg);
        fprintf('\t\tPercent of pixels you get to keep = %.4f (# = %i)\n',1-pd/(n1+n2),pd)
    end
end


if other_indices_flag
    fprintf('\tComputing other distance measures...\n')

    %GT_L = bwlabel(GT_L,4);
    %Alg_L = bwlabel(Alg_L,4);

    L1s = reshape(L_gt,num_rows*num_cols,1) + 1;
    L2s = reshape(L_alg,num_rows*num_cols,1) + 1;
    [AR,RI,MI,HI] = valid_RandIndex(L1s,L2s);
    fprintf('\t\tRand index: %.3f\n',RI)
    fprintf('\t\tAdjusted Rand index: %.3f\n',AR)
    fprintf('\t\tMirkin index: %.3f\n',MI)
    fprintf('\t\tHubert index: %.3f\n',HI)
end

toc;

end



function [L_mistakes,correct,incorrect] = boundary_dist(L1,L2)
%BOUNDARY_DIST computes the distance between boundaries.
%
% Measures how many boundary edges in L1 are preserved in L2.
%
% Algorithm: for each boundary pixel p in L1, look at the neighborhood
% around p in L2 and see if it contains a boundary pixel. If so, mark p has
% matched/correct. If not, mark it as incorrect.
%
% Because the ground truth borders do not align with the watershed borders
% exactly, we define a loose boundary around p to which we match p to.
%
% Returns:
%   L_mistakes = label matrix of incorrect pixels.
%   correct    = # of boundary pixels correctly matched.
%   incorrect  = # of boundary pixels not matched/incorrect.
%


% size of the neighborhood/loose border.
nhood_radius = 2;

[num_rows num_cols] = size(L1);
L_mistakes = zeros(size(L1)); % matrix of 1's where mistakes occurred.
[boundary_rows, boundary_cols] = find(L1==0);
correct = 0;
incorrect = 0;


for i=1:size(boundary_rows,1)

    % get boundary pixel coordinate.
    x = boundary_rows(i);
    y = boundary_cols(i);

    % get loose border around the boundary pixel, but in L2.
    t = L2(max(1,x-nhood_radius):min(num_rows,x+nhood_radius), max(1,y-nhood_radius):min(num_cols,y+nhood_radius));

    % check if it's near a boundary.
    if ismember(0,t) == 1
        correct = correct + 1;
    else
        incorrect = incorrect + 1;
        L_mistakes(x,y) = 1;
    end

end


end



function M = compute_m(L1,L2,connsize)
%COMPUTE_M computes the matching matrix, M.
%
% M(i,j) is the intersection of region i from L1 and region j from L2.
%

[L1,num_regions1] = bwlabel(L1,4);
[L2,num_regions2] = bwlabel(L2,connsize);

Region1 = regionprops(L1,'PixelIdxList');
Region2 = regionprops(L2,'PixelIdxList');

M = zeros(num_regions1,num_regions2);

fprintf('\t\t# of regions: L1=%i, L2=%i\n',num_regions1, num_regions2)


%% compute the M matrix.
for i=1:num_regions1
    for j=1:num_regions2
        M(i,j) = size(intersect_sorted(Region1(i).PixelIdxList, Region2(j).PixelIdxList),1);
        %M(i,j) = size(intersect(Region1(i).PixelIdxList, Region2(j).PixelIdxList),1);
    end
end


end


function [n2,asym_dist] = asym_dist(M,L1,L2)
% ASYM_DISTANCE computes the asymmetric partition distance between two
% partitions.
%
% The asymmetric partition distance is the minimum number of elements that
% must be deleted from L2 so that L2 is finer than L1. Finer means that
% every region in L2 is a subset of exactly one region in L1. So this
% measures gives a sense of how many boundary pixels are crossed.
%
% Algorithm: n2 - \sum_i (max_j M(i,j))
%   where n2 = the number of non-zero pixels in L2.
%         M(i,j) = |intersect(L1_i,L2_j)|
%
% From: Toward a Generic Evaluation of Image Segmentation
%       Jaime Cardoso and Luis Corte-Real
%       Transactions on Image Processing, Nov 2005.
%
% Example:
%
%   L1 = 1 1 0 0 2      L2 = 1 1 0 0 2      L3 = 1 1 1 1 1
%        0 0 0 2 2           1 0 0 2 2           1 0 0 1 1
%        3 3 0 0 0           0 0 0 0 0           0 0 0 0 0
%        3 3 0 4 0           3 3 0 0 4           2 2 0 0 3
%        3 3 0 4 4           3 3 0 4 4           2 2 0 3 3
%
%  >> gt_asym_dist(L1,L2) = 2
%  >> gt_asym_dist(L2,L1) = 3
%  >>
%
% In our case, L1 is the ground truth and L2 is the auto segmentation.
%

n2 = size(find(L2>0),1); % # of non-zero pixels.
asym_dist = n2 - sum(max(M));

end



function [n1,n2,partition_dist] = partition_dist(M,L1,L2)
% PARTITION_DIST computes the partition distance between the two
% partitions.
%
% The partition distance is the number of pixels that must be deleted from
% L1 and L2 such that L1 exactly matching L2. If the number of non-zero
% pixels is the same in L1 and L2 (a typical partition on the same data),
% then the partition distance is: n - matching_cost, where n is the number
% of elements in the dataset.
%
% In our case, however, the number of non-zero pixels can be different in
% L1 and L2. Instead, we return: n1 + n2 - 2*matching_cost, where n1 and n2
% are the number of non-zero pixels in L1 and L2, respectively.
%
% From: Partition-Distance: A Problem and Class of Perfect Graphs Arising
%                           in Clustering
%       Dan Gusfield
%       Information Processing Letters, May 2002
%
% Example:
%
%   L1 = 1 1 0 0 2      L2 = 1 1 0 0 2      L3 = 1 1 1 1 1
%        0 0 0 2 2           1 0 0 2 2           1 0 0 1 1
%        3 3 0 0 0           0 0 0 0 0           0 0 0 0 0
%        3 3 0 4 0           3 3 0 0 4           2 2 0 0 3
%        3 3 0 4 4           3 3 0 4 4           2 2 0 3 3
%
%  >> partition_dist(L1,L2) = 5
%  >> partition_dist(L1,L3) = 11
%
% This distance is symmetric; the lower the better.
% Interpretation: partition_dist / (n1+n2) is the percentage of total nodes
% that have to be deleted.
%


% we want min-matching, so transform the problem.
max_val = max(max(M));
Mt = abs(M-max_val);


%% solve the matching problem in Mt.
[rowsol,cost,v,u,costMat] = lapjv(Mt);

% rowsol contains the indices of the matched partners. go through and
% comput the cost wrt M.
cost = 0;
for i=1:size(rowsol,2)
    if rowsol(i) == 0
        continue
    else
        % rowsol is wrt the larger of the two assignments halves.
        if size(Mt,1) > size(Mt,2)
            cost = cost + M(rowsol(i),i);
        else
            cost = cost + M(i,rowsol(i));
        end
    end
end

% compute the partition distance.
n1 = size(find(L1>0),1);
n2 = size(find(L2>0),1);
partition_dist = n1 + n2 - 2*cost;

end



function boundary_dist2(L1,L2)

[rows1 cols1] = find(L1==0);
[rows2 cols2] = find(L2==0);

[num_rows num_cols] = size(L2);

n1 = size(rows1,1)
n2 = size(rows2,1)

% give a label to each boundary pixel in L1 and L2.
L1_label = zeros(size(L1));
l = 1;
for i=1:n1
    x = rows1(i);
    y = cols1(i);

    L1_label(x,y) = l;
    l = l + 1;
end

L2_label = zeros(size(L2));

l = 1;
for i=1:n2
    x = rows2(i);
    y = cols2(i);

    L2_label(x,y) = l;
    l = l + 1;
end

toc;

M = sparse(n1,n2);
d_max = 10;

% if there exists a boundary pixel within dist, then add that to the match.

nhood_radius = 5;

% go through boundary pixels p in L1.
for i=1:n1
    x = rows1(i);
    y = cols1(i);

    % look in search window around p in L2.
    for r=max(1,x-nhood_radius):min(num_rows,x+nhood_radius)
        for s=max(1,y-nhood_radius):min(num_cols,y+nhood_radius)
            if L2_label(r,s) ~= 0
                d = sqrt((x-r)^2 + (y-s)^2);
                if d < d_max
                    M(L1_label(x,y),L2_label(r,s)) = d;
                    %M(L2_label(r,s),L1_label(x,y)) = d;
                end
            end
        end
    end
end

toc;

fprintf('Computing matching\n')
%[rowsol,cost,v,u,costMat] = lapjv(M);
[matching,cost] = munkres(M);
cost

toc;

end


%=========================================================================%
%=========================================================================%

function validation_teshting()
%% teshting.
L1 = [1 1 0 0 2; 0 0 0 2 2; 3 3 0 0 0; 3 3 0 4 0; 3 3 0 4 4];
L2 = [1 1 0 0 2; 1 0 0 2 2; 0 0 0 0 0; 3 3 0 0 4; 3 3 0 4 4];
L3 = [1 1 1 1 1; 1 0 0 1 1; 0 0 0 0 0; 3 3 0 0 4; 3 3 0 4 4];


M12 = compute_m(L1,L2);
M21 = compute_m(L2,L1);
M13 = compute_m(L1,L3);
M31 = compute_m(L3,L1);
M23 = compute_m(L2,L3);
M32 = compute_m(L3,L2);


%% partition_dist.
x = partition_dist(M12,L1,L2);
if x ~= 5, error('Partition distance between L1,L2 is %i not %i',x,5), end

x = partition_dist(M21,L2,L1);
if x ~= 5, error('Partition distance between L2,L1 is %i not %i',x,5), end

x = partition_dist(M13,L1,L3);
if x ~= 11, error('Partition distance between L1,L3 is %i not %i',x,11), end

x = partition_dist(M31,L3,L1);
if x ~= 11, error('Partition distance between L3,L1 is %i not %i',x,11), end

x = partition_dist(M23,L2,L3);
if x ~= 8, error('Partition distance between L2,L3 is %i not %i',x,8), end

x = partition_dist(M32,L3,L2);
if x ~= 8, error('Partition distance between L3,L2 is %i not %i',x,8), end


%% asym_dist.
x = asym_dist(M12,L1,L2);
if x ~= 2, error('Asym distance between L1 and L2 is %i not %i',x,2), end

x = asym_dist(M21,L2,L1);
if x ~= 3, error('Asym distance between L2 and L1 is %i not %i',x,3), end

x = asym_dist(M13,L1,L3);
if x ~= 6, error('Asym distance between L1 and L3 is %i not %i',x,6), end

x = asym_dist(M31,L3,L1);
if x ~= 3, error('Asym distance between L3 and L1 is %i not %i',x,3), end

x = asym_dist(M32,L3,L2);
if x ~= 0, error('Asym distance between L3 and L2 is %i not %i',x,0), end

x = asym_dist(M23,L2,L3);
if x ~= 5, error('Asym distance between L2 and L3 is %i not %i',x,5), end



end




% GRAVEYARD
% munkres.
%toc;
%[matching,cost] = munkres(Mt);
%cost = 0;
%for i=1:size(matching,2)
%    if matching(i) == 0
%        continue
%    else
%        cost = cost + M(i,matching(i));
%    end
%end
%toc;
%cost_alg1 = cost;

%hungarian.
%[matching,cost] = Hungarian(Mt);
% compute the cost wrt the original intersection size numbers (the
% non-transformed matrix).
%cost = sum(M(matching==1));
%toc;
