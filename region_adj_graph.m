%========================================================================%
%    Region adjacency graph                                              %
%                                                                        %
%    Author: Saket Navlakha                                              %
%    Date  : August 2010                                                 %
%                                                                        %
%========================================================================%


function [Edges, Boundary_size, Boundary] = region_adj_graph(L)
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

Boundary_size = sparse(num_regions,num_regions); % [reg1,reg2] -> # of shared boundary pixels.
Boundary = struct(); % e.g. Boundary(reg).Pixels = [x1 y1; x2 y2; ...]

% allocate space for the boundary pixels.
Regions = regionprops(L,'Area');
for i=1:num_regions
    Boundary(i).Pixels = zeros(3*Regions(i).Area,2);
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
            Boundary_size(reg1,reg2) = Boundary_size(reg1,reg2) + 1;
            Boundary_size(reg2,reg1) = Boundary_size(reg2,reg1) + 1;
            
            % store the boundary pixels for each region.
            t = curr_boundary_row(reg1);
            Boundary(reg1).Pixels(t,:) = [x y];
            curr_boundary_row(reg1) = t + 1;
            
            t = curr_boundary_row(reg2);
            Boundary(reg2).Pixels(t,:) = [x y];
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
            Boundary_size(reg1,reg2) = Boundary_size(reg1,reg2) + 1;
            Boundary_size(reg2,reg1) = Boundary_size(reg2,reg1) + 1;
            
            % store the boundary pixels for each region.
            t = curr_boundary_row(reg1);
            Boundary(reg1).Pixels(t,:) = [x y];
            curr_boundary_row(reg1) = t + 1;
            
            t = curr_boundary_row(reg2);
            Boundary(reg2).Pixels(t,:) = [x y];
            curr_boundary_row(reg2) = t + 1;            
                       
        % 8hood cases: don't store edges; only boundaries. schematic:
        %   1 2
        %    p
        %   3 4 
        % TODO: just look at all remaining 4 directions at once and add
        % them pairwise.
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
                Boundary_size(reg1,reg2) = Boundary_size(reg1,reg2) + 1;
                Boundary_size(reg2,reg1) = Boundary_size(reg2,reg1) + 1;
                 
                % store the boundary pixels for each region.
                t = curr_boundary_row(reg1);
                Boundary(reg1).Pixels(t,:) = [x y];
                curr_boundary_row(reg1) = t + 1;
            
                t = curr_boundary_row(reg2);
                Boundary(reg2).Pixels(t,:) = [x y];
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
    Boundary(i).Pixels(start_row:end_row,:) = [];
end


end