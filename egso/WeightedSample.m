function idx = WeightedSample(w)
% WeightedSample  Draw one index according to the weight vector w (no Statistics Toolbox needed)
%
% Equivalent to: randsample(length(w), 1, true, w)
% i.e. draws one index with probability proportional to w
%
% Input:   w    positive weight vector (normalised internally)
% Output:  idx  integer index, 1 <= idx <= length(w)

w   = w(:)' ;           % row vector
w   = w / sum(w);       % normalise to probabilities
cdf = cumsum(w);        % cumulative distribution
r   = rand();           % uniform random number in [0,1)
idx = find(cdf >= r, 1, 'first');   % inverse-transform sampling
if isempty(idx)
    idx = length(w);    % defensive fallback
end
end
