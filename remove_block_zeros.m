%========================================================================%
%    REMOVES BLOCKS OF ZEROS                                             %
%                                                                        %
%    Author: Saket Navlakha                                              %
%    Date  : August 2010                                                 %
%                                                                        %
%========================================================================%

function L = remove_block_zeros(L)
%REMOVE_BLOCK_ZEROS     Removes squares of zeros.
%
% REMOVE_BLOCK_ZEROS(L) removes square blocks of zeros that do not 
% themselves separate two regions.
%
% For example, this:
%
%   1 1 1 1
%   1 0 0 1
%   1 0 0 1
%   1 1 1 1
%
% ...becomes:
%
%   1 1 1 1
%   1 1 1 1
%   1 1 1 1
%   1 1 1 1
%


[num_rows, num_cols] = size(L);


for x=1:num_rows,
    for y=1:num_cols,

        % ignore non-boundary pixels.
        if L(x,y) ~= 0
            continue
        end

        if (x+2 > num_rows || y+2 > num_cols || x-1 < 1 || y-1 < 1)
            continue
        end

        % at north west (top left), check if east, south, and south-east are 0.
        if (L(x,y+1) == 0 && L(x+1,y) == 0 && L(x+1,y+1) == 0)

            % check if all pixels surrounding the block are the same.
            w = L(x-1,y); % start north.

            if (w == 0)
                continue
            elseif (L(x-1,y+1) == w && L(x-1,y+2) == w && L(x,y+2) == w && L(x+1,y+2) == w && L(x+2,y+2) == w && L(x+2,y+1) == w && L(x+2,y) == w && L(x+2,y-1) == w && L(x+1,y-1) == w && L(x,y-1) == w && L(x-1,y-1) == w)

                % convert the square to w.
                L(x,y) = w;
                L(x,y+1) = w;
                L(x+1,y) = w;
                L(x+1,y+1) = w;
            end

        end
    end
end

    

end