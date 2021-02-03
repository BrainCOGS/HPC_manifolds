function [BafterA_r, BafterA_l, BafterA_r_shf, BafterA_l_shf] = shuffleDoublets_v2( nshuffles, alltrials_choice, lefttrials_choice, righttrials_choice, temp3)

ntrials = length(alltrials_choice);

ncells = size(temp3{1},2);

A = zeros(2, ncells, ntrials); %will just put together max and min times into a 3D array, and work with that

t = temp3;
for itrial = 1:ntrials
    
    a = t{itrial};
    nt  = size(a,1); %may want to change this to a vector
    
    sa = sum(a);
    fi = find(sa ~=0);
    
    an = a(:,fi); %non zero spiking cells only
    
    %but better to convert to times
    at = bsxfun(@times, an, (1:nt)');
    
    at(at==0) = nan; %probably bad practice but should get rid of the zeros so that the min doesn't give zeros
    
    % On what trials does B spike after A?
    % I.e. on each trial, is there an element for which max(B) > min(A)
    
    maxb = max(at);
    mina = min(at);
  
    A(:, fi, itrial) = [maxb; mina];     
end


%now everything should be straightforward - can shuffle the array A, not
%worrying about the different trial lengths, different numbers of spikes

Al = A(:,:, lefttrials_choice);
Ar = A(:,:, righttrials_choice);
nltrials = length(lefttrials_choice);
nrtrials = length(righttrials_choice);


BafterA_l_shf = zeros(ncells,ncells,nshuffles);
BafterA_r_shf = zeros(ncells,ncells,nshuffles);

for ishf = 1:nshuffles

    
    Als = zeros(size(Al));
    Ars = zeros(size(Ar));

    %for each shuffle we can shuffle the trials per cell.
    
    tordl = zeros(ncells, nltrials); %track the new order, just in case
    for ic = 1:ncells
        tordl(ic,:) = randperm(nltrials);
        Als(:,ic,:) = Al(:,ic, tordl(ic,:)); %rearrange trials for cell ic only
    end
    
    tordr = zeros(ncells, nrtrials); %track the new order, just in case
    for ic = 1:ncells
        tordr(ic,:) = randperm(nrtrials);
        Ars(:,ic,:) = Ar(:,ic, tordr(ic,:)); %rearrange trials for cell ic only
    end
    
    Alsmin = permute(Als(2,:,:),[2 1 3]);
    Alsmax = Als(1,:,:);
    
    Arsmin = permute(Ars(2,:,:),[2 1 3]);
    Arsmax = Ars(1,:,:);
    
    bls = bsxfun(@and, bsxfun( @minus, Alsmax, Alsmin) > 0, Alsmin>0); %make sure Als had a spike
        
    BafterA_l_shf(:,:,ishf) = sum(bls,3);
       
    brs = bsxfun(@and, bsxfun( @minus, Arsmax, Arsmin) > 0, Arsmin>0); %make sure Ars had a spike
    
    BafterA_r_shf(:,:,ishf) = sum(brs,3);
    
    disp(['Shuffle # ' num2str(ishf) ])
    
end

%real data
Almin = permute(Al(2,:,:),[2 1 3]);
Almax = Al(1,:,:);

Armin = permute(Ar(2,:,:),[2 1 3]);
Armax = Ar(1,:,:);

% bl = bsxfun( @minus, Almax, Almin) > 0;
bl = bsxfun(@and, bsxfun( @minus, Almax, Almin) > 0, Almin>0);

BafterA_l = sum(bl,3);

% br = bsxfun( @minus, Armax, Armin) > 0;
br = bsxfun(@and, bsxfun( @minus, Armax, Armin) > 0, Armin>0);

BafterA_r = sum(br,3);




end





