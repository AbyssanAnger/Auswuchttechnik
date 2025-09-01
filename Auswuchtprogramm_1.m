load("awt_UE1_Aufgabe.mat");
load("awt_UE1_Zusatz1.mat");
mean_value = mean(awt_UE1_Zusatz1(:, 4))*10^-6;
test_freq = 1/mean_value;
% figure;
% for i = 1:size(awt_UE1, 2)
%     subplot(size(awt_UE1, 2), 1, i);
%     plot(awt_UE1(:, i));
%     title(['Column ' num2str(i)]);
%     xlabel('Index');
%     ylabel(['Value of Column ' num2str(i)]);
% end
% Identifying the rising edge of the 3rd column in awt_UE1
column_data = awt_UE1_Zusatz1(:, 3);
threshold = mean(column_data); % Define a threshold for rising edge detection
rising_edges = find(diff(column_data > threshold) == 1); % Find indices of rising edges

% Cutting the 1st, 2nd, and 3rd columns at the first and last rising edge
cut_indices = [rising_edges(1)+1, rising_edges(end)];
cut_data = awt_UE1_Zusatz1(cut_indices(1):cut_indices(2), 1:3);

% Performing FFT for each column in cut_data
fft_results = zeros(size(cut_data)); % Preallocate for FFT results
frequencies = cell(1, 3); % Cell array to store frequency vectors

for k = 1:3
    fft_results(:, k) = fft(cut_data(:, k))/length(cut_data); % Compute FFT
    frequencies{k} = linspace(0, 1, length(cut_data)) * (1 / (1 / size(cut_data, 1))); % Frequency vector
end

% Generating a single-sided amplitude spectrum from the FFT results
single_sided_amplitude_spectrum = cell(1, 3); % Cell array to store single-sided amplitude spectra

for k = 1:3
    L = length(cut_data); % Length of the cut data
    P2 = abs(fft_results(:, k)); % Two-sided spectrum
    P1 = P2(1:L/2+1); % Single-sided spectrum
    P1(2:end-1) = 2*P1(2:end-1); % Multiply by 2 for all but the DC and Nyquist components
    single_sided_amplitude_spectrum{k} = P1; % Store the single-sided amplitude spectrum
end

single_sided_amplitude_spectrum_array = cell2mat(single_sided_amplitude_spectrum);

% Halving the frequency range for the single-sided amplitude spectrum
for k = 1:3
    half_index = 1:floor(length(single_sided_amplitude_spectrum{k})/2); % Indices for halving
    single_sided_amplitude_spectrum{k} = single_sided_amplitude_spectrum{k}(half_index); % Halve the amplitude spectrum
    halfed_frequencies{k} = frequencies{k}(half_index); % Halve the frequency vector
end

% Finding the highest frequency after 0 in the single-sided amplitude spectrum of the 3rd column
[~, max_index] = max(single_sided_amplitude_spectrum{3}); % Find index of maximum amplitude
highest_frequency_after_zero = halfed_frequencies{3}(max_index); % Get the corresponding frequency
disp(['The highest frequency after 0 in the SSAS of the 3rd column is: ' num2str(highest_frequency_after_zero) ' Hz']);

% Finding the FFT values from the first two columns at the position of max_index
fft_values_at_max_index = zeros(2, 1); % Preallocate for FFT values
for j = 1:2
    fft_values_at_max_index(j) = single_sided_amplitude_spectrum_array(max_index, j); % Get FFT values at max_index for columns 1 and 2
end
disp(['FFT values at the max index are: ' num2str(fft_values_at_max_index')]);

% Finding the FFT angles from the first two columns at the position of max_index
fft_angles_at_max_index = zeros(2, 1); % Preallocate for FFT values
for j = 1:2
    fft_angles_at_max_index(j) = angle(fft_results(max_index, j))*180/pi; % Get FFT values at max_index for columns 1 and 2
end
disp(['FFT angle at the max index columns are: ' num2str(fft_angles_at_max_index')]);

% Plotting the single-sided amplitude spectrum (SSAS) for the 1st, 2nd, and 3rd columns
figure;
for n = 1:3
    subplot(3, 1, n);
    plot(halfed_frequencies{n}, single_sided_amplitude_spectrum{n});
    title(['Single-Sided Amplitude Spectrum - Column ' num2str(n)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
end

% Creating a polar plot with FFT angles and FFT values
figure;
for m = 1:2
    subplot(2, 1, m);
    polarplot(deg2rad(fft_angles_at_max_index(m)), fft_values_at_max_index(m), 'o');
    title(['Polar Plot - Column ' num2str(m)]);
    thetalim([0 360]); % Set theta limits
    rlim([0 max(fft_values_at_max_index)]); % Set radius limits
end

% % Plotting the FFT results for the 1st, 2nd, and 3rd columns
% figure;
% for m = 1:3
%     subplot(3, 1, m);
%     plot(frequencies{m}, abs(fft_results(:, m)));
%     title(['FFT Results - Column ' num2str(m)]);
%     xlabel('Frequency (Hz)');
%     ylabel('Magnitude');
% 
% end

% % Plotting the FFT results with frequency, amplitude, and phase angle in a 3D plot
% figure;
% for m = 1:3
%     subplot(3, 1, m);
%     % Create a 3D plot for each column
%     plot3(frequencies{m}, abs(fft_results(:, m)), angle(fft_results(:, m)));
%     title(['3D FFT Results - Column ' num2str(m)]);
%     xlabel('Frequency (Hz)');
%     ylabel('Magnitude');
%     zlabel('Phase Angle (radians)');
%     grid on;
% end

% Plotting the cut data for the 1st, 2nd, and 3rd columns
figure;
for j = 1:3
    subplot(3, 1, j);
    plot(cut_data(:, j));
    title(['Cut Data - Column ' num2str(j)]);
    xlabel('Index');
    ylabel(['Value of Column ' num2str(j)]);
end


% Preallocate STFT results
stft_results = cell(1, 3); % Cell array to store STFT results

for k = 1:3
    [S, F, T] = stft(cut_data(:, k),test_freq,Window=kaiser(1500),OverlapLength=1400,FFTLength=1500); % Compute STFT
    stft_results{k} = S; % Store the STFT results for the current column
end

% Plotting the spectrogram for each column
figure;
for k = 1:3
    subplot(3, 1, k);
    mesh(T, F, abs(stft_results{k})); % Plot spectrogram
    title(['Spectrogram - Column ' num2str(k)]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    view(2);
    colorbar;
end

figure;
mesh(T, F, abs(stft_results{1})./F); % Plot spectrogram
    title(['Spectrogram - Column ' num2str(3)]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    view(2);
    ylim([0 50]);
    clim([0 10]);
    colorbar;