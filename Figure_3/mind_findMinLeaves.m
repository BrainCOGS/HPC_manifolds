function outLeaves = mind_findMinLeaves(dat,xInputs)
    % transforms distance data in the MIND dataset
    % such that powerlaw can be fitted.
    
    % First: rebin and take only the relevant i->j distances
    Dtrial = dat.forestdat.rwd.Dg;
    upperTriangularMatrix = boolean(triu(ones(size(Dtrial)),1));
    distances = Dtrial(upperTriangularMatrix);
    
    % Second, count the relevant section
    Cs = zeros(length(xInputs),1);
    for idx = 1:length(xInputs)
       Cs(idx) = sum( distances<xInputs(idx) );
    end

    % Third: save data, and some metadata
    outLeaves.Cs = Cs;
    outLeaves.Dtrial = dat.forestdat.rwd.Dg;
    outLeaves.Nframes = length(dat.forestdat.data.f);

end