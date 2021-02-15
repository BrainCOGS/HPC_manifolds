function normOut = norm_columns(matIn)
    
    % Outputs a row vector whose components are the L2 norm of
    % the columns of matIn
    
    assert(numel(size(matIn)) == 2); % only accepts 2d arrays (ie, matrices)
    normOut = (sum(matIn.^2,1)).^(1/2);
    
end