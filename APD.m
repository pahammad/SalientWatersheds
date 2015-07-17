%========================================================================%
%    ASYMMETRIC PARTITION DISTANCE                                       %
%                                                                        %
%    Author: Saket Navlakha                                              %
%    Date  : June 7, 2011                                                %
%                                                                        %
%========================================================================%


function dist = APD(L1,L2)
% APD computes the asymmetric partition distance between two partitions.
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
%  >> L1 = [1 1 0 0 2; 0 0 0 2 2; 3 3 0 0 0; 3 3 0 4 0; 3 3 0 4 4];
%  >> L2 = [1 1 0 0 2; 1 0 0 2 2; 0 0 0 0 0; 3 3 0 0 4; 3 3 0 4 4];
%  >> L3 = [1 1 1 1 1; 1 0 0 1 1; 0 0 0 0 0; 3 3 0 0 4; 3 3 0 4 4];
%
%   L1 = 1 1 0 0 2      L2 = 1 1 0 0 2      L3 = 1 1 1 1 1
%        0 0 0 2 2           1 0 0 2 2           1 0 0 1 1
%        3 3 0 0 0           0 0 0 0 0           0 0 0 0 0
%        3 3 0 4 0           3 3 0 0 4           2 2 0 0 3
%        3 3 0 4 4           3 3 0 4 4           2 2 0 3 3
%
%  >> APD(L1,L2) = 2
%  >> APD(L2,L1) = 3
%  >> APD(L1,L3) = 6
%  >> APD(L3,L1) = 3
%  >> APD(L2,L3) = 5
%  >> APD(L3,L2) = 0
%
% In our case, L1 is the ground truth and L2 is the algorithm's partition.
%


%fprintf('\tComputing M...\n')
M = compute_m(L1,L2);

%fprintf('\tComputing APD...\n')
n = size(find(L2>0),1); % # of non-zero pixels.
dist = (n - sum(max(M))) / n; %uncomment for percentage

end



function M = compute_m(L1,L2)
%COMPUTE_M computes the matching matrix, M.
%
% M(i,j) is the intersection (number of common pixels) of region i from L1 
% and region j from L2.
%
% TODO: can ignore first steps if we assume bwlabel is done beforehand.

[L1,num_regions1] = bwlabel(L1,4);
[L2,num_regions2] = bwlabel(L2,4);

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



function c = intersect_sorted(a,b)
%INTERSECT_SORTED     Set intersection between sorted sets.
% INTERSECT_SORTED(A,B) when A and B are vectors returns the values common
% to both A and B.  A and B must be sorted and unique, and the result will be
% sorted and unique.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

c = a(ismember_sorted(a,b));

end



function [tf,loc] = ismember_sorted(a,s)
%ISMEMBER_SORTED   True for member of sorted set.
% ISMEMBER_SORTED(A,S) for the vector A returns an array of the same size as A
% containing 1 where the elements of A are in the set S and 0 otherwise.
% A and S must be sorted and cannot contain NaN.
%
% [TF,LOC] = ISMEMBER_SORTED(A,S) also returns an index array LOC where
% LOC(i) is the index in S which matches A(i) (highest if there are ties)
% or 0 if there is no such index.
%
% See also ISMEMBER, MATCH_SORTED, INTERSECT_SORTED, SETDIFF_SORTED, UNION_SORTED.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

% The internal function ismembc comes from ismember.m
% It requires non-sparse arrays.
a = full(a);
s = full(s);
if nargout < 2
  tf = ismembc(a,s);
else
  loc = ismembc2(a,s);
  tf = (loc > 0);
end

end

end