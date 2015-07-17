%========================================================================%
%    Post-processing routines                                            %
%                                                                        %
%    Author: Saket Navlakha                                              %
%    Date  : August 2010                                                 %
%                                                                        %
%========================================================================%

function Q = post_process(Q,xx)
%POST_PROCESS post processes the label matrix.
%
% POST_PROCESS(Q,xx) removes small regions of size < 'xx' and leftover
% boundary regions (such as squares).
%

%xx = 3;
Q = remove_boundaries(Q);
Q = remove_small_regions(Q,xx);
Q = remove_boundaries(Q);
Q = remove_block_zeros(Q);

end
