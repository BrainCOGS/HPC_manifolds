function nci = nieh_ci(idata, numBoot)

% Makes the 95% CI for shaded error bars

for i=1:size(idata,2)
    
    idata1 = idata(:,i);
    
    nci(:,i) = bootci(numBoot,{@mean,idata1});
    
end