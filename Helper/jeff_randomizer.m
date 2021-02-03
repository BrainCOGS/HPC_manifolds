function dFF_random = jeff_randomizer(dataDFF, jeff_constant)

% Changed T = sum(window) to make sure shift can't be too small (7/31/2020)
% set RNG (7/31/2020)

rng(42);

[T,N] = size(dataDFF);
all = 0:(jeff_constant/T):jeff_constant;

for t = 1:jeff_constant
    window = ((t-1) <= all) & (all < t);
    T = sum(window);
    for i = 1:N
        Tshift = floor(T/4)+randi(floor(T/2));
        dFF_random(window,i) = circshift( dataDFF(window,i), Tshift);
    end
end



end