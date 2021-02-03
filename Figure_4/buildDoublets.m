function [BafterA, ball] = buildDoublets(temp3)

t = cellfun(@(x) logical(x), temp3, 'uniformoutput',false);

ntrials = length(temp3);
nc = size(temp3{1},2);

ball = false(nc,nc, ntrials);

for itrial = 1:ntrials
    
    a = t{itrial};
    nt  = size(a,1); %may want to change this to a vector
    
    sa = sum(a);
    fi = find(sa ~=0);
    
    an = a(:,fi); %non zero spiking cells only
    
    %but better to convert to times
    at = bsxfun(@times, an, (1:nt)');
    
    at(at==0) = nan; %probably bad practice but should get rid of the zeros so that the min doesn't give zeros
    
    
    % I understood you to be calculating:
    % On what trials does B spike after A?
    % I.e. on each trial, is there an element for which max(B) > min(A)
    
    maxb = max(at);
    mina = min(at);
    
    b = bsxfun( @minus, maxb, mina') > 0; %if you want spikes in the same bin to be ok then change this to >=
    
    ball(fi, fi, itrial) = b;
end

%%
% if you want a list for each cell pair, then I would do this at the end
% but this is slower than the whole previous stage, because of the for loop

BafterA = cell(nc,nc);
for ii = 1:nc
    for ij = 1:nc
        BafterA{ii,ij} = find(ball(ii,ij,:));
    end
end








