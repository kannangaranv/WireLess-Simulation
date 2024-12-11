M = 4; % Number of transmit antennas
N = 4; % Number of receive antennas
Pt_vals = [1, 2, 5, 10]; 
sigma2 = 1; % Noise variance
num_iterations = 1000; 

rates_equal = zeros(length(Pt_vals), 1);
rates_inversion = zeros(length(Pt_vals), 1);
rates_best_eigenmode = zeros(length(Pt_vals), 1);
rates_waterfilling = zeros(length(Pt_vals), 1);

for pt_idx = 1:length(Pt_vals)
    Pt = Pt_vals(pt_idx); 
    rate_eq = 0;
    rate_inv = 0;
    rate_best = 0;
    rate_wf = 0; 
    
    for iter = 1:num_iterations
        % Generate random channel matrix H
        H = (randn(N, M) + 1j * randn(N, M)) / sqrt(2); 
        
        % SVD decomposition
        [~, S, ~] = svd(H);
        eigen_values = diag(S).^2; 
        
        %----------------------------------------------------------------
        % Equal Power Allocation
        p_eq = Pt / M; 
        rate_eq = rate_eq + sum(log2(1 + (p_eq * eigen_values) / sigma2));
        
        %-----------------------------------------------------------------
        % Channel Inversion Power Allocation
        p_inv = ((1 ./ eigen_values) / sum(1 ./ eigen_values)) * Pt;
        rate_inv = rate_inv + sum(log2(1 + p_inv .* eigen_values / sigma2));
        
        %-----------------------------------------------------------------
        % Best Eigenmode Allocation
        rate_best = rate_best + log2(1 + Pt * max(eigen_values) / sigma2);
        
        %-----------------------------------------------------------------
        % Waterfilling Power Allocation 
        alpha = eigen_values / sigma2; 
        
        M_wf = M; % Number of active streams in waterfilling
        Pk = zeros(1, M);

        while M_wf > 0
            lambda_inv = (sum(1 ./ alpha(1:M_wf)) + Pt) / M_wf; % 1/lambda
            Pk(1:M_wf) = lambda_inv - 1 ./ alpha(1:M_wf); % Power allocation
            
            % Check for negative power allocation
            if all(Pk(1:M_wf) >= 0)
                break; 
            else
                [~, minIdx] = min(Pk(1:M_wf));
                alpha(minIdx) = [];
                M_wf = M_wf - 1; 
            end
        end
        
        % Compute achievable rate 
        rate = 0;
        for i = 1:M_wf
            rate = rate + (log2(1 + Pk(i) .* alpha(i) / sigma2)); 
        end
        
        rate_wf = rate_wf + rate; 
    end
    
    % Average rates over all iterations
    rates_equal(pt_idx) = rate_eq / num_iterations;
    rates_inversion(pt_idx) = rate_inv / num_iterations;
    rates_best_eigenmode(pt_idx) = rate_best / num_iterations;
    rates_waterfilling(pt_idx) = rate_wf / num_iterations;

    % Display the average rates for each scheme at the current transmit power Pt
    disp(['Transmit Power: ', num2str(Pt), ' W']);
    disp(['Equal Power Allocation Rate: ', num2str(rates_equal(pt_idx)), ' bps/Hz']);
    disp(['Channel Inversion Allocation Rate: ', num2str(rates_inversion(pt_idx)), ' bps/Hz']);
    disp(['Best Eigenmode Allocation Rate: ', num2str(rates_best_eigenmode(pt_idx)), ' bps/Hz']);
    disp(['Waterfilling Allocation Rate: ', num2str(rates_waterfilling(pt_idx)), ' bps/Hz']);
    disp('----------------------------------------------------');

end

% Plot results
figure;
plot(Pt_vals, rates_equal, '-o', 'DisplayName', 'Equal Power');
hold on;
plot(Pt_vals, rates_inversion, '-s', 'DisplayName', 'Channel Inversion');
plot(Pt_vals, rates_best_eigenmode, '-^', 'DisplayName', 'Best Eigenmode');
plot(Pt_vals, rates_waterfilling, '-x', 'DisplayName', 'Waterfilling');
grid on;
xlabel('Transmit Power (W)');
ylabel('Achievable Rate (bps/Hz)');
title('Achievable Rate for Different Power Allocation Schemes');
legend('Location', 'Best');
