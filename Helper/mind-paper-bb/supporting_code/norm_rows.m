function normOut = norm_rows(matIn)
    
    % Outputs a row vector whose components are the L2 norm of
    % the rows of matIn
    normOut = transpose( norm_columns( transpose( matIn ) ) );
    
end