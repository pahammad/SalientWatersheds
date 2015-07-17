%========================================================================%
%    Agglomerative region merging                                        %
%                                                                        %
%    Author: Saket Navlakha                                              %
%                                                                        %
%    v0.2, 05/24/2011 -- added #SPs stopping criteria                    %
%    v0.3, 06/09/2011 -- modified merging criteria to use weighted       %
%                        canny edges                                     %
%                                                                        %
%========================================================================%

function Lx = AgglomerateBSPs(I, L, S, stopSPs)
%AGGLOMERATEBSPS    Runs the region merging algorithm.
%
% AGGLOMERATEBSPS(I_nlm,I,L,S,numSPs) runs the region merging algorithm given
% I (the original image), L (the label matrix), S (the weighted salient edges),
% and numSPs (the final number of SPs).
%
% Returns: a label matrix with the specified number of regions at several
% points in the hierarchical tree decomposition.
%


%% Flags, constants, and paths.

tic;

addpath('texture/');
addpath('FastEMD/');

boundary_size_flag = 0; % TODO: leftover from old code.
intensity_flag = 1; % use intensity features.
texture_flag = 1;   % use texture features.
salient_flag = 0;   % use salient edges. [[replaced boundary_size_flag.]]
alpha = 1;          % weight of the texture component.

num_histogram_bins = 32;    % # of bins in the intensity/texture histograms.
emd_threshold = 8;  % EMD distance threshold.
EMD_dist = ones(num_histogram_bins+1,num_histogram_bins+1).*emd_threshold;
for i=1:num_histogram_bins+1
    for j=max([1 i-emd_threshold+1]):min([num_histogram_bins+1 i+emd_threshold-1])
        EMD_dist(i,j)= abs(i-j); 
    end
end

[num_rows num_cols] = size(L);


%% Get adjacent regions.
fprintf('\t\tGetting adjacent regions and region boundaries...\n')
[Edges, Shared_Boundary_Size, Boundary_Pixels] = region_adj_graph(L);
num_regions = max(max(L));

fprintf('\t\t# of regions: %i\n',num_regions)
fprintf('\t\t# of edges  : %i\n',size(Edges,1))



%% Compute pixel properties.
fprintf('\t\tComputing pixel properties...\n')

Regions = regionprops(L,I,'PixelValues','PixelIdxList'); % for intensity and texture.

if (intensity_flag)
    intensity_bins = 1:(256/num_histogram_bins):257;
    Reg2Intensity = zeros(num_regions,num_histogram_bins+1);  % region -> pixel intensity histogram.
end

if (texture_flag) 
    I_texture = MR8fast(I)';
    
    % get max and min values along each dimension. this is used to make
    % equally-sized and spaced-histograms.
    max_texture = max(I_texture);
    min_texture = min(I_texture);
    
    if size(max_texture,2) ~= 8 || size(min_texture,2) ~= 8
        error('Expected 8 texture dimensions but found %i,%i',size(max_texture,2),size(min_texture,2))
    end
    
    texture_bins = zeros(num_histogram_bins+1,8);
    for i=1:8
        texture_bins(:,i) = (min_texture(i) : (max_texture(i)-min_texture(i))/num_histogram_bins : max_texture(i));
    end

    % for each region: 8 histograms, each of size 'num_histogram_bins+1'.
    Reg2Texture = zeros(num_regions,num_histogram_bins+1,8);
end


if (intensity_flag && texture_flag)
    if size(intensity_bins,2) ~= size(texture_bins,1)
        error('Size of intensity bins (%ix%i) does not equal size of texture bins (%ix%i)',size(intensity_bins,1),size(intensity_bins,2),size(texture_bins,1),size(texture_bins,2))
    end
end


%% Grouping pixel-level properties to the region level.
fprintf('\t\tComputing region properties...\n')
for reg=1:num_regions
    
    if (intensity_flag)
        h = histc(Regions(reg).PixelValues,intensity_bins); % edge centers.
        %h = h./sum(h); % convert to probability distribution.
        % 02/01/2013: replaced above line with the if-else below for the
        % case where there's only 0-intensity pixels in the region. The
        % above line caused a Nan.
        if sum(h) == 0
            h = h;
        else
            h = h./sum(h); % convert to probability distribution.
        end
        Reg2Intensity(reg,:) = h;

        if isnan(h), error('Intensity histogram for region %i is NaN',reg), end        
        if sum(h) ~= 0 && abs(sum(h)-1) > 0.00001, error('Intensity histogram does not sum to 1: %i',reg), end
    end
    
    if (texture_flag)        
        region_texture = I_texture(Regions(reg).PixelIdxList,:); % # region pixels x 8.
        
        for i=1:8
            h = histc(region_texture(:,i),texture_bins(:,i));
            h = h ./ sum(h);
            Reg2Texture(reg,:,i) = h;
            
            if isnan(h), error('Texture histogram for region %i is NaN',reg), end
            if abs(sum(h)-1) > 0.00001, error('Texture histogram does not sum to 1: %i',reg), end
        end        
    end
end


%% Create salient region adjacency matrix, Sx.

if (salient_flag)
    fprintf('\tComputing shared salient boundaries...\n')
    
    S = S ./ max(S(:)); % normalize.
    
    % go thru S and compute the sum of the salient weights on the boundary between two regions.
    Shared_Salient_Weights = sparse(num_regions,num_regions); % reg1,reg2 -> sum of salient weights on shared boundary.
    
    for x=1:num_rows,
        for y=1:num_cols
            
            % salient edge and watershed boundary.
            if (S(x,y) > 0 && L(x,y) == 0)
                
                % north and south.
                if (x-1 > 0) && (x+1 <= num_rows) && (L(x-1,y) ~= 0) && (L(x+1,y) ~= 0) && (L(x-1,y) ~= L(x+1,y))
                    Shared_Salient_Weights(L(x-1,y),L(x+1,y)) = Shared_Salient_Weights(L(x-1,y),L(x+1,y)) + S(x,y);
                    Shared_Salient_Weights(L(x+1,y),L(x-1,y)) = Shared_Salient_Weights(L(x+1,y),L(x-1,y)) + S(x,y);
                end
                
                % west and east.
                if (y-1 > 0) && (y+1 <= num_cols) && (L(x,y-1) ~= 0) && (L(x,y+1) ~= 0) && (L(x,y-1) ~= L(x,y+1))
                    Shared_Salient_Weights(L(x,y-1),L(x,y+1)) = Shared_Salient_Weights(L(x,y-1),L(x,y+1)) + S(x,y);
                    Shared_Salient_Weights(L(x,y+1),L(x,y-1)) = Shared_Salient_Weights(L(x,y+1),L(x,y-1)) + S(x,y);
                end
                
                % 8hood.
                if (x-1 > 0 && x+1 <= num_rows) && (y-1 > 0 && y+1 <= num_cols)
                    % northeast and southwest.
                    if L(x-1,y+1) ~= 0 && L(x+1,y-1) ~= 0 && L(x-1,y+1) ~= L(x+1,y-1)
                        Shared_Salient_Weights(L(x-1,y+1),L(x+1,y-1)) = Shared_Salient_Weights(L(x-1,y+1),L(x+1,y-1)) + S(x,y);
                        Shared_Salient_Weights(L(x+1,y-1),L(x-1,y+1)) = Shared_Salient_Weights(L(x+1,y-1),L(x-1,y+1)) + S(x,y);
                    end

                    % northwest and southeast.
                    if L(x+1,y+1) ~= 0 && L(x-1,y-1) ~= 0 && L(x+1,y+1) ~= L(x-1,y-1)
                        Shared_Salient_Weights(L(x+1,y+1),L(x-1,y-1)) = Shared_Salient_Weights(L(x+1,y+1),L(x-1,y-1)) + S(x,y);
                        Shared_Salient_Weights(L(x-1,y-1),L(x+1,y+1)) = Shared_Salient_Weights(L(x-1,y-1),L(x+1,y+1)) + S(x,y);
                    end
                    
                end
            end            
        end
    end

    clear S;

end


%% Create heap matrix of pairwise costs/weighted adj region graph, H.
fprintf('\t\tComputing H...\n')

% set variables to use.
if (~intensity_flag), Reg2Intensity = []; end
if (~texture_flag), Reg2Texture = []; end
if (~salient_flag), Shared_Boundary_Size = []; Shared_Salient_Weights = []; end

num_regions = double(num_regions); % added for matlab 2011.

H = sparse(num_regions,num_regions);

for i=1:size(Edges)
    reg1 = Edges(i,1);
    reg2 = Edges(i,2);
    
    if (reg1 == reg2), error('Adj matrix contains self-loop: %i, %i', reg1,reg2), end

    F1 = Reg2Intensity(reg1,:);
    F2 = Reg2Intensity(reg2,:); 
    
    % The higher the value of d, the closer the pair (because we do the exp thing).    
    d = region_distance(reg1,reg2,alpha,EMD_dist,Reg2Intensity,Reg2Texture,Shared_Salient_Weights,Shared_Boundary_Size);
    %d = region_distance(reg1,reg2,EMD_dist,Reg2Intensity,Reg2Texture,Shared_Boundary_Size,Regions);
       
    if (d < 0 || d > 1 || isnan(d)), error('Distance not within bound for %i,%i: %.2f', reg1,reg2,d), end

    % Standard approach, no boosting of small regions.
    %H(reg1,reg2) = d;
    %H(reg2,reg1) = d;

    % Added 10/26/2012 to boost small regions so they are merged faster.
    min_size_boost = exp(-min(length(Regions(reg1).PixelIdxList),length(Regions(reg2).PixelIdxList))/10);    
        
    H(reg1,reg2) = d+min_size_boost;
    H(reg2,reg1) = d+min_size_boost;
    
end


clear Edges;
%fprintf('\t\t# of heap entries: %i\n',size(find(H>0),1)/2)


%% Merge adjacent regions if their merger abides by the opt function.
fprintf('\t\tMerging...\n')

while (1)

    % go through heap and find max cost pair.
    [colmax,colmax_ind] = max(H);
    [rowcolmax,rowcolmax_ind] = max(colmax);
    
    best_cost = rowcolmax(1);
    reg1 = rowcolmax_ind;
    reg2 = colmax_ind(rowcolmax_ind);
       
    % shtatus update and save.
    if mod(num_regions,500) == 0
        fprintf('\t\t\t@ %i,%.4f\n',num_regions,str2double(num2str(best_cost)))
        
        % save the output.
        if num_regions <= 12000
            Lx.(strcat('m',int2str(num_regions))) = L; 
        end
        
    end
        
    if num_regions == stopSPs
        fprintf('\t\tReached desired number of regions: %i. Exiting...\n',num_regions)
        Lx.(strcat('m',int2str(num_regions))) = L;
        break;
    end
    
    %if (numSPs == num_regions)
    %    fprintf('\t\tReached desired number of regions: %i. Exiting...\n',num_regions)
    %    ASPIm = L;
    %    break;
    
    % error checking.
    %if (best_cost > 1 || isnan(best_cost))
    if (best_cost > 2 || isnan(best_cost)) % with the min_size_boost, the best_cost could be greater than 1.
        Lx.(strcat('m',int2str(num_regions))) = L;
        error('Invalid best cost: %.3f for regions %i,%i. Exiting...\n',str2double(num2str(best_cost)),reg1,reg2)
    elseif (reg1 == reg2)
        Lx.(strcat('m',int2str(num_regions))) = L;
        error('Trying to merge the same region: %i, %i',reg1,reg2)
    end
       
    % merge regions, remove boundary pixels.
    old_label = max(reg1,reg2);
    new_label = min(reg1,reg2);
    [L,Boundary_Pixels,Regions] = merge_regions(L,Boundary_Pixels,Regions,old_label,new_label);

    
    % update new_label region properties, delete old_label.
    new_label_pixels = Regions(new_label).PixelIdxList;    
    if (intensity_flag)
        Reg2Intensity(old_label,:) = 0;
        h = histc(I(new_label_pixels),intensity_bins);
        %h = histc(I_nlm(new_label_pixels),intensity_bins);
        h = h./sum(h);
        Reg2Intensity(new_label,:) = h;
        
        if isnan(h), error('Intensity histogram for region %i is NaN',new_label), end
        if abs(sum(h)-1) > 0.00001, error('Intensity histogram does not sum to 1: %i',new_label), end        
    end

    if (texture_flag)        
        Reg2Texture(old_label,:,:) = 0;
        region_texture = I_texture(new_label_pixels,:);
        for i=1:8
            h = histc(region_texture(:,i),texture_bins(:,i));
            h = h ./ sum(h);
            Reg2Texture(new_label,:,i) = h;  
            %Reg2Texture(reg,:,i) = h;
            
            if isnan(h), error('Texture histogram for region %i is NaN',new_label), end
            if abs(sum(h)-1) > 0.00001, error('Texture histogram does not sum to 1: %i',new_label), end            
        end        
    end

    
    % update the heap and the new salient edge boundaries.       
    old_label_neighbors = find(H(:,old_label)>0); % 'find' returns indices.
    for n = old_label_neighbors'

        if (n == old_label), error('Self-loop detected: %i, %i',n,old_label), end

        % if the salient edge exists between the old_label and the
        % neighbor, then the salient edge now exists between new_label and
        % the neighbor. (e.g. before Sx(3,5) = 0 but if Sx(4,5) = 1 then
        % change Sx(3,5) = 1.
        if (new_label ~= n)
            if (salient_flag && Sx(new_label,n) == 0)            
                if (Sx(old_label,n) ~= Sx(n,old_label)), error('Non-symmetry found in Sx: %i,%i',old_label,n), end

                Sx(new_label,n) = Sx(old_label,n);
                Sx(n,new_label) = Sx(old_label,n);
            end
        end
        
        % update the graph, part i: add old_label to the list of
        % new_label's neighbors and vice-versa.
        if (new_label ~= n)
            if H(old_label,n) ~= H(n,old_label), error('Non-symmetry found in H: %i,%i',old_label,n), end
                        
            % just temporarily to ensure they are nonzero.
            H(new_label,n) = 100;
            H(n,new_label) = 100;
        end
        
        % update the boundary sizes.
        % TODO: are there corner cases here? probably, but ignore for now.
        if (new_label ~= n)
            if (boundary_size_flag)
                if (Shared_Boundary_Size(old_label,n) <= 0)
                    error('Boundary size between two neighboring regions is 0: %i,%i,%i',old_label,n,Boundary_size(old_label,n))
                end
                Shared_Boundary_Size(new_label,n) = Shared_Boundary_Size(new_label,n) + Shared_Boundary_Size(old_label,n);
                Shared_Boundary_Size(n,new_label) = Shared_Boundary_Size(n,new_label) + Shared_Boundary_Size(old_label,n);
            end
        end
        
    end

    % update the graph, part ii: remove old_label from the graph.
    H(:,old_label) = 0;
    H(old_label,:) = 0;
    if (boundary_size_flag)
        Shared_Boundary_Size(:,old_label) = 0;
        Shared_Boundary_Size(old_label,:) = 0;
    end
    num_regions = num_regions - 1;

    % update the heap: for all neighbors of new_label, compute their cost.
    new_label_neighbors = find(H(:,new_label)>0);
    for n = new_label_neighbors'
        
        % check if the pair lie between a salient edge.
        if (salient_flag && Sx(new_label,n) == 1)
            % small and non-zero so neighbors can be found.
            H(new_label,n) = 0.001;
            H(n,new_label) = 0.001;
        else          
            d = region_distance(new_label,n,alpha,EMD_dist,Reg2Intensity,Reg2Texture,Shared_Salient_Weights,Shared_Boundary_Size);
            %d = reg_distance(new_label,n,EMD_dist,Reg2Intensity,Reg2Texture,Shared_Boundary_Size,Regions);
            
            if (d < 0 || d > 1 || isnan(d))
                error('Distance not within bound for %i,%i: %.2f', new_label,n,d)
            end
            
            

        % Standard approach, no boosting of small regions.
        %H(new_label,n) = d;
        %H(n,new_label) = d;

        % Added 10/26/2012 to boost small regions so they are merged faster.
        min_size_boost = exp(-min(length(Regions(new_label).PixelIdxList),length(Regions(n).PixelIdxList))/10);
        
        H(new_label,n) = d+min_size_boost;
        H(n,new_label) = d+min_size_boost;
    
        end
    end
    
end

toc;

end




function [Edges, Shared_Boundary_Size, Boundary_Pixels] = region_adj_graph(L)
%REGION_ADJ_GRAPH returns the edges in the region adj graph defined by 'L'.
%
% REGION_ADJ_GRAPH(L) computes the region adjacency matrix of label image
% 'L'. Two regions are adjacent if they are separated by a 0 pixel in the
% horizontal or vertical direction.
%
% Returns: 
%   -Edges: two-column matrix of region adjacencies sorted in ascending
%   order.
%   -Boundary_size: a (num_regions,num_regions)-sized matrix with entry
%   (i,j) = the number of pixels defining the boundary between i and j.
%   -Boundary: a structure mapping each region to its boundary pixels.
%       ( e.g. Boundary(reg).Pixels = [x1 y1; x2 y2; ...] )
%
% Edges are defined with 4hood. Boundary* are defined with 8hood.
%
% Example:
%
% img =
%
%     1     1     1     0     3
%     0     0     0     0     0
%     4     0     0     2     0
%     4     4     0     2     2
%     0     0     0     0     0
%     5     5     5     5     5
%
% >> [Edges,Boundary_size, Boundary] = region_adj_graph(img)
%
% Edges =
%
%   1   3
%   1   4
%   2   4
%   2   5
%   4   5
%

[num_rows num_cols] = size(L);
num_regions = max(max(L));

% TODO: what's the maximum number of edges in a planar graph?
Edges = zeros(1000000,2); % reg1, reg2
curr_edge_row = 1;

num_regions=double(num_regions);
Shared_Boundary_Size = sparse(num_regions,num_regions); % [reg1,reg2] -> # of shared boundary pixels.
Boundary_Pixels = struct(); % e.g. Boundary(reg).Pixels = [x1 y1; x2 y2; ...]

% allocate space for the boundary pixels.
Regions = regionprops(L,'Area');
for i=1:num_regions
    Boundary_Pixels(i).Pixels = zeros(3*Regions(i).Area,2);
end
curr_boundary_row = ones(num_regions,1); % current boundary row for each region.


%% find region adjacencies.
for x=1:num_rows,
    for y=1:num_cols,
        
        if L(x,y) ~= 0
            continue
        end
        
        % north and south.
        if (x-1 > 0) && (x+1 <= num_rows) && (L(x-1,y) ~= 0) && (L(x+1,y) ~= 0) && (L(x-1,y) ~= L(x+1,y))
            
            reg1 = L(x-1,y);
            reg2 = L(x+1,y);
            
            % populate the edges.
            Edges(curr_edge_row,:) = [reg1 reg2];
            curr_edge_row = curr_edge_row + 1;
            
            % increment boundary size between the two regions.
            Shared_Boundary_Size(reg1,reg2) = Shared_Boundary_Size(reg1,reg2) + 1;
            Shared_Boundary_Size(reg2,reg1) = Shared_Boundary_Size(reg2,reg1) + 1;
            
            % store the boundary pixels for each region.
            t = curr_boundary_row(reg1);
            Boundary_Pixels(reg1).Pixels(t,:) = [x y];
            curr_boundary_row(reg1) = t + 1;
            
            t = curr_boundary_row(reg2);
            Boundary_Pixels(reg2).Pixels(t,:) = [x y];
            curr_boundary_row(reg2) = t + 1;
            

        % North and west (or south and east) can't be nonzero and different
        % because of "consecutive_nonzeros.m". Hence, the "elseif".        
        % west and east.
        elseif (y-1 > 0) && (y+1 <= num_cols) && (L(x,y-1) ~= 0) && (L(x,y+1) ~= 0) && (L(x,y-1) ~= L(x,y+1))
            
            reg1 = L(x,y-1);
            reg2 = L(x,y+1);
            
            % populate the edges.
            Edges(curr_edge_row,:) = [reg1 reg2];
            curr_edge_row =  curr_edge_row + 1;
            
            % increment boundary size between the two regions.
            Shared_Boundary_Size(reg1,reg2) = Shared_Boundary_Size(reg1,reg2) + 1;
            Shared_Boundary_Size(reg2,reg1) = Shared_Boundary_Size(reg2,reg1) + 1;
            
            % store the boundary pixels for each region.
            t = curr_boundary_row(reg1);
            Boundary_Pixels(reg1).Pixels(t,:) = [x y];
            curr_boundary_row(reg1) = t + 1;
            
            t = curr_boundary_row(reg2);
            Boundary_Pixels(reg2).Pixels(t,:) = [x y];
            curr_boundary_row(reg2) = t + 1;            
        
        % 8hood cases: don't store edges; only boundaries. schematic:
        %   1 2
        %    p
        %   3 4 
        % TODO: just look at all remaining 4 directions at once and add
        % them pairwise.
        % TODO: 06/09/2011 -- not sure right now why we need this part.
        % leaving it for now.
        elseif (y-1 > 0 && y+1 <= num_cols) && (x-1 > 0 && x+1 <= num_rows)
            
            to_add = zeros(6,2); % [reg1 reg2]
            add_row = 1;
            
            % 1 and 2.
            if L(x-1,y-1) ~= 0 && L(x-1,y+1) ~= 0 && L(x-1,y-1) ~= L(x-1,y+1)
                to_add(add_row,:) = [L(x-1,y-1), L(x-1,y+1)];
                add_row = add_row + 1;
            end
            
            % 1 and 3.
            if L(x-1,y-1) ~= 0 && L(x+1,y-1) ~= 0 && L(x-1,y-1) ~= L(x+1,y-1)
                to_add(add_row,:) = [L(x-1,y-1), L(x+1,y-1)];
                add_row = add_row + 1;
            end
            
            % 1 and 4.
            if L(x-1,y-1) ~= 0 && L(x+1,y+1) ~= 0 && L(x-1,y-1) ~= L(x+1,y+1)
                to_add(add_row,:) = [L(x-1,y-1), L(x+1,y+1)];
                add_row = add_row + 1;
            end

            % 2 and 3.
            if L(x-1,y+1) ~= 0 && L(x+1,y-1) ~= 0 && L(x-1,y+1) ~= L(x+1,y-1)
                to_add(add_row,:) = [L(x-1,y+1), L(x+1,y-1)];
                add_row = add_row + 1;
            end
            
            % 2 and 4.
            if L(x-1,y+1) ~= 0 && L(x+1,y+1) ~= 0 && L(x-1,y+1) ~= L(x+1,y+1)
                to_add(add_row,:) = [L(x-1,y+1), L(x+1,y+1)];
                add_row = add_row + 1;
            end

            % 3 and 4.
            if L(x+1,y-1) ~= 0 && L(x+1,y+1) ~= 0 && L(x+1,y-1) ~= L(x+1,y+1)
                to_add(add_row,:) = [L(x+1,y-1), L(x+1,y+1)];
                add_row = add_row + 1;
            end
            
            % uniquify and add to Boundary variables.
            to_add(add_row:6,:) = []; % delete unused rows.
            to_add = unique(sort(to_add,2),'rows');
            
            for i=1:size(to_add,1)
                reg1 = to_add(i,1);
                reg2 = to_add(i,2);
                
                %increment boundary size between the two regions.
                Shared_Boundary_Size(reg1,reg2) = Shared_Boundary_Size(reg1,reg2) + 1;
                Shared_Boundary_Size(reg2,reg1) = Shared_Boundary_Size(reg2,reg1) + 1;
                 
                % store the boundary pixels for each region.
                t = curr_boundary_row(reg1);
                Boundary_Pixels(reg1).Pixels(t,:) = [x y];
                curr_boundary_row(reg1) = t + 1;
            
                t = curr_boundary_row(reg2);
                Boundary_Pixels(reg2).Pixels(t,:) = [x y];
                curr_boundary_row(reg2) = t + 1;     
            end            
        end      
        
        if curr_edge_row > size(Edges,1)
            error('Not enough space preallocated for adjacency edges.\n')
        end
    end
end


%% sort edges by increasing value, remove duplicates.

% 'sort' ensures that n1<n2, and increasing order of n2 for n1=constant.
% necessary to do this first because otherwise [1 2] will be treated
% separately from [2 1].
%edges = sortrows(sort(edges, 2));
Edges = sort(Edges, 2);

% remove duplicate edges
Edges = unique(Edges, 'rows');

% remove the first row because it has 0 0 from the preallocation.
% TODO: does this require recopying the entire matrix?
Edges(1,:) = [];


%% remove zeros from the Boundary pixels.
for i=1:num_regions
    start_row = curr_boundary_row(i);
    end_row = 3*Regions(i).Area;
    Boundary_Pixels(i).Pixels(start_row:end_row,:) = [];
end


end



function dist = region_distance(reg1,reg2,alpha,EMD_dist,Reg2Intensity,Reg2Texture,Shared_Salient_Weights,Shared_Boundary_Size)
%REGION_DISTANCE The merging distance between two regions.
%
% dist = foo*exp(EMD(Int(reg1),Int(reg2)) +  \alpha * \sum_{i=1}^8 EMD(Text(reg1,i),Text(reg2,i)) / 8)),
%       where foo = sum of salient weights along shared boundary / # shared boundary pixels
%       where alpha = the weight of the texture component.
%
% Then, the function returns: exp(-d)
%
% To not use intensity, set Reg2Intensity = []. Similar for texture
% (Reg2Texture) and salient edges (Shared_Salient_Weights).
%

dist = 0;

if size(Reg2Intensity,1) ~= 0
    
    F1 = Reg2Intensity(reg1,:);
    F2 = Reg2Intensity(reg2,:);    
    dist = emd_hat_gd_metric_mex(F1',F2',EMD_dist);
end

if size(Reg2Texture,1) ~= 0
    
    F1 = Reg2Texture(reg1,:,:);
    F2 = Reg2Texture(reg2,:,:);    
    t = 0;
    
    for i=1:8
        F1i = F1(1,:,i);
        F2i = F2(1,:,i);
        
        t = t + emd_hat_gd_metric_mex(F1i',F2i',EMD_dist);
    end
    
    dist = dist + (alpha*t/8);
end

% TODO: how to use foo??

dist = exp(-dist);


if size(Shared_Salient_Weights,1) ~= 0
    if Shared_Salient_Weights(reg1,reg2) ~= Shared_Salient_Weights(reg2,reg1)
        error('Non-symmetry in Shared_Salient_Weights for regions %i,%i: %i,%i',reg1,reg2,Shared_Salient_Weights(reg1,reg2),Shared_Salient_Weights(reg2,reg1));
    end
    
    if Shared_Boundary_Size(reg1,reg2) ~= Shared_Boundary_Size(reg2,reg1)
        error('Non-symmetry in Shared_Boundary_Size for regions: %i,%i: %i,%i',reg1,reg,Shared_Boundary_Size(reg1,reg2),Shared_Boundary_Size(reg2,reg1));
    end
    
    foo = Shared_Salient_Weights(reg1,reg2)/Shared_Boundary_Size(reg1,reg2)+1e-10; % just non-zero so the other terms don't vanish.
    
    if (foo < 0 || foo > 1 || isnan(foo))
        error('foo issues: %s',foo); 
    end
    
    dist = foo * dist;   
end 

end