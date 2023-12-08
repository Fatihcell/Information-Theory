function entropy = DiffEnt(data)
        
    % This function gives differential entropy based on gaussian kernel PDF estimation. 
    % data: input data vector
    
    % Created by Fatih Onay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n         = length(data); 
    pdf       = zeros(1, m); % Assign empty PDF 
    npnts     = 0.1*m; % Number of states
    xi        = linspace(min(data), max(data), npnts); % Points to evaluate KDE
    bandwidth = 1.06 * std(x) * length(x)^(-1/5); % bandwidth selection
    m         = length(xi);

    for i = 1:m
        % Calculate the kernel values for each data point at xi(i)
        kernels = exp(-0.5 * ((xi(i) - data) / bandwidth).^2) / (bandwidth * sqrt(2 * pi));
        
        % Sum and normalize the kernels to get the PDF
        pdf(i) = sum(kernels) / n;
    end
     entropy = -sum(pdf .* log2(pdf + eps)) * mean(diff(xi));
end
