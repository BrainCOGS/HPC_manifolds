function [landmarks] = get_landmarks(data,n)
    % Returns the indices of 'n' landmarks in the dataset 'data' using a 
    % MaxMin algorithm with greedy optimization.
    %
    % Landmark points are chosen one at a time, and each new landmark 
    % maximizes, over all unused data points, the minimum  euclidian 
    % distances between states to any of the existing landmarks. 
    % The first point is chosen arbitrarily. [De Silva & Tenenbaum 2004].

    data = data(1:end-1,:); % all states are viable landmarks except 
                            % the last point because it has no followers
                            
    [T,N] = size(data);
    landmarks = zeros(n,1); 
    m = zeros(T,1); % contains, for every state, the distance to the next landmark point.
    landmarks(1) = randi(T,1,1);
        
    for j = 1:T
        if j ~= landmarks(1)
            m(j) = sum( (data(landmarks(1),:) - data(j,:)).^2 );
        else
            m(j) = nan;
        end
    end

    for i = 2:n
        [v, ind] = max(m);
        landmarks(i) = ind;
        for j = 1:T
            s = sum( (data(ind,:)-data(j,:)).^2  );
            if s == 0
                m(j) = nan; % this point is already used.
            elseif m(j) > s
                m(j) = s;   % update m only, if j is closer to the new landmark than all others.
            end
        end
        %disp([num2str(100*i/n) '% of landmarks found.'])
    end
    landmarks = sort(landmarks);
end

