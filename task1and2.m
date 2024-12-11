Eb_N0_dB = -10:1:10;            
Eb_N0 = 10.^(Eb_N0_dB / 10);  
num_bits = 10^4;  
N = 1000;

%------------------------------------------------------------------
% Task 1: Average bit error rate (BER) evaluation in fading channels

%------------------------------------------------------------------
% a)Simulation:

sim_ber_rayleigh = zeros(size(Eb_N0_dB)); 
sim_ber_noise = zeros(size(Eb_N0_dB));

for ii = 1:length(Eb_N0)
    BER_fading=0;
    BER_noise=0;
    % Generate random BPSK symbols
    for jj = 1:N
        bits = randi([0, 1], 1, num_bits);
        symbols = sqrt(Eb_N0(ii))*(2 * bits - 1); 

        % Generate Rayleigh fading channel
        hh = sqrt(1/2) * (randn(1, num_bits) + 1j * randn(1, num_bits)); 

        % Add AWGN noise
        noise = sqrt(1/2) * (randn(1, num_bits) + 1j * randn(1, num_bits));
        received = symbols.*hh + noise; 
        received_noise = symbols + noise;
        
        % Equalize the channel
        equalized = conj(hh).*received;
        
        % Detect bits
        detectedBits = equalized > 0;
        detectedNoise = received_noise > 0;
        BER_fading = BER_fading + mean(abs(detectedBits-bits));
        BER_noise = BER_noise + mean(abs(detectedNoise-bits));
    end
    % Compute BER
    sim_ber_rayleigh(ii) = BER_fading/N;
    sim_ber_noise(ii) = BER_noise/N;
end

% --------------------------------------------------------------
% c) Theoretical BER for Rayleigh Fading:

theoretical_ber_rayleigh = (1 - sqrt(Eb_N0 ./ (1 + Eb_N0))) / 2;

% --------------------------------------------------------------
% d) Theoretical BER for AWGN

theoretical_ber_awgn = qfunc(sqrt(2 * Eb_N0)); 

% --------------------------------------------------------------
% Task 2: Comparing diversity combining schemes

M = 3; 
BER_sim_MRC = zeros(length(Eb_N0_dB), 1);
BER_theory_MRC = zeros(length(Eb_N0_dB), 1);
BER_sim_SC = zeros(length(Eb_N0_dB), 1);

for ii = 1:length(Eb_N0_dB)
    BER_fading_MRC = 0;
    %------------------------------------------------------------------
    %  a) Simulation of MRC receiver

    % Generate random BPSK symbols
    for jj = 1:N
        bits = randi([0, 1], 1, num_bits);
        symbols = sqrt(Eb_N0(ii))*(2*bits - 1);

        % Generate Rayleigh fading channels for each branch
        hh = sqrt(1/2)*(randn(M, num_bits) + 1i*randn(M, num_bits));

        % Add noise to each branch (AWGN)
        noise = sqrt(1/2)*(randn(M, num_bits) + 1i*randn(M, num_bits));

        receivedSignal = symbols.*hh + noise;
        
        %-------------------------------------------------
        % a) MRC: Sum the received signals from all branches
        combinedSignal = sum(conj(hh) .* receivedSignal, 1);

        % Demodulate
        detectedBits = combinedSignal > 0;

        BER_fading_MRC = BER_fading_MRC + mean(abs(detectedBits-bits));
       
    
    end
    % Calculate the Bit Error Rate (BER)
    BER_sim_MRC(ii) = BER_fading_MRC / N;
    
    % -----------------------------------------------------
    % b) Theoretical BER calculation 
    
    gamma = Eb_N0(ii); 
    Gamma = sqrt(gamma / (1 + gamma));
    summation = 0;

    for m = 0:M-1
        summation = summation + nchoosek(M-1+m, m) * ( (1 + Gamma)/2 )^m;
    end
    
    BER_theory_MRC(ii) = (((1 - Gamma)/2)^M) * summation;
    
end

%-----------------------------------------------------------------------
% c) Selection Combining: Pick the branch with the maximum signal-to-noise ratio (SNR)

for ii = 1:length(Eb_N0_dB)
    BER_fading = 0; 

    for jj = 1:N
        
        bits = randi([0, 1], 1, num_bits);

        symbols = sqrt(Eb_N0(ii)) * (2 * bits - 1);

        hh = sqrt(1/2) * (randn(M, num_bits) + 1i * randn(M, num_bits));
        
        noise = sqrt(1/2) * (randn(M, num_bits) + 1i * randn(M, num_bits));

        receivedSignal = symbols .* hh + noise;

        % Find the branch with the maximum SNR (highest |h|^2)
        [~, bestBranch] = max(abs(hh).^2, [], 1);

        % Select the signal and channel corresponding to the best branch
        selectedSignal = receivedSignal(sub2ind(size(receivedSignal), bestBranch, 1:num_bits));
        selectedChannel = hh(sub2ind(size(hh), bestBranch, 1:num_bits));
        
        equalized = conj(selectedChannel) .* selectedSignal;

        detectedBits = real(equalized) > 0;

        BER_fading = BER_fading + mean(detectedBits ~= bits);
    end

    BER_sim_SC(ii) = BER_fading / N;
end

% Plot Results
figure;
semilogy(Eb_N0_dB, sim_ber_rayleigh, 'b-o', 'LineWidth', 1.5); hold on;
semilogy(Eb_N0_dB, theoretical_ber_rayleigh, 'g-', 'LineWidth', 1.5);
semilogy(Eb_N0_dB, sim_ber_noise, 'm-o', 'LineWidth', 1.5);
semilogy(Eb_N0_dB, theoretical_ber_awgn, 'r-', 'LineWidth', 1.5);
semilogy(Eb_N0_dB, BER_sim_MRC, 'k-o', 'LineWidth', 1.5);
semilogy(Eb_N0_dB, BER_theory_MRC, 'y-', 'LineWidth', 1.5);
semilogy(Eb_N0_dB, BER_sim_SC, 'm-', 'LineWidth', 1.5);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER of BPSK in Rayleigh Fading');
legend('Simulated Rayleigh', 'Theoretical Rayleigh', 'Theoretical AWGN','Simulated Noise', 'Simulated BER of MRC','Theoretical BER of MRC ', 'Simulated BER of SC');