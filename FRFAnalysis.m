
%% FRF analysis
%Setting of parameters
%fs= 4000;% TODO CHECK OF WAT OP LINE 48 STaat hiervoor ook werkt. Zodat dit puur parameters zijn voor anaylise.
nfft=16000%resolutie = fs/nfft
noverlap = nfft/2;
window= hann(nfft);

% Reducing the size of the parameters to test
%y=y(1:length(y)/2);
%u=u(1:length(u)/2);

% Getting the estimate
[H,hz] = tfestimate(u, y, window, noverlap, nfft, fs);

% Making the magnitude/phase plot
H_mag = abs(H);
H_db = mag2db(H_mag);
H_angle_rad = angle(H);
H_angle_deg = rad2deg(H_angle_rad);

subplot(2,1,1)
semilogx(hz,H_db)
title("Bode plot of FRF")
ylabel("Gain (dB)")

subplot(2,1,2) 
semilogx(hz,H_angle_deg)
%title("Phase H")
ylabel("Phase (Degrees)")
xlabel("Frequency (s^{-1})")
%xlabel("Angular frequency (s^{-1})")

%subplot(3,1,3)
%[cxy, w]=mscohere(u, y, window, noverlap, nfft, fs);
%semilogx(w, cxy)
%grid on
%hold on
%title("Coherence")
%xlabel("Frequency (s^{-1})")
%ylabel("Coherence")

%loglog(hz, abs(H))

%% Coherance plots
nfft=8000%resolutie = fs/nfft
noverlap = nfft/2;
window= hann(nfft);
[cxy, w]=mscohere(u, y, window, noverlap, nfft, fs);
semilogx(w, cxy, 'LineWidth', 2)
%grid on
hold on
title("Coherence")
xlabel("Frequency (s^{-1})")
ylabel("Coherence")


%% Aditional coherance aditions nfft
subplot(2,1,1)
nfft=40000;%resolutie = fs/nfft
noverlap = nfft/2;
window= hann(nfft);
%subplot(3,1,3)
[cxy, w]=mscohere(u, y, window, noverlap, nfft, fs);
semilogx(w, cxy)
%grid on
hold on
title("Coherence for different resolutions")
xlabel("Frequency (s^{-1})")
ylabel("Coherence")

subplot(2,1,1)
nfft=16000%resolutie = fs/nfft
noverlap = nfft/2;
window= hann(nfft);
[cxy, w]=mscohere(u, y, window, noverlap, nfft, fs);
semilogx(w, cxy, 'LineWidth', 2)

subplot(2,1,1)
nfft=8000;%resolutie = fs/nfft
noverlap = nfft/2;
window= hann(nfft);
%subplot(3,1,3)
[cxy, w]=mscohere(u, y, window, noverlap, nfft, fs);
semilogx(w, cxy)

subplot(2,1,1)
nfft=4000%resolutie = fs/nfft
noverlap = nfft/2;
window= hann(nfft);
[cxy, w]=mscohere(u, y, window, noverlap, nfft, fs);
semilogx(w, cxy)

subplot(2,1,1)
nfft=2000;%resolutie = fs/nfft
noverlap = nfft/2;
window= hann(nfft);
%subplot(3,1,3)
[cxy, w]=mscohere(u, y, window, noverlap, nfft, fs);
semilogx(w, cxy)

subplot(2,1,1)
nfft=1000;%resolutie = fs/nfft
noverlap = nfft/2;
window= hann(nfft);
%subplot(3,1,3)
[cxy, w]=mscohere(u, y, window, noverlap, nfft, fs);
semilogx(w, cxy)

legend("0.1 Hz", "0.25 HZ", "0.5 HZ", "1 Hz", "2 Hz", "4 Hz")%, "1/64 Data")

%% Aditional coherance aditions amount of data
subplot(2,1,2)
nfft=16000%resolutie = fs/nfft
noverlap = nfft/2;
window= hann(nfft);
[cxy, w]=mscohere(u, y, window, noverlap, nfft, fs);
semilogx(w, cxy, 'LineWidth', 2)
%grid on
hold on
title("Coherence for fractions of data")
xlabel("Frequency (s^{-1})")
ylabel("Coherence")

subplot(2,1,2)
y_copy=y(1:length(y)/2);
u_copy=u(1:length(u)/2);
%subplot(3,1,3)
[cxy, w]=mscohere(u_copy, y_copy, window, noverlap, nfft, fs);
semilogx(w, cxy)

subplot(2,1,2)
y_copy=y(1:length(y)/4);
u_copy=u(1:length(u)/4);
%subplot(3,1,3)
[cxy, w]=mscohere(u_copy, y_copy, window, noverlap, nfft, fs);
semilogx(w, cxy)

subplot(2,1,2)
y_copy=y(1:length(y)/8);
u_copy=u(1:length(u)/8);
%subplot(3,1,3)
[cxy, w]=mscohere(u_copy, y_copy, window, noverlap, nfft, fs);
semilogx(w, cxy)

subplot(2,1,2)
y_copy=y(1:length(y)/64);
u_copy=u(1:length(u)/64);
%subplot(3,1,3)
%[cxy, w]=mscohere(u_copy, y_copy, window, noverlap, nfft, fs);
%semilogx(w, cxy)

legend("All data", "1/2 Data", "1/4 Data", "1/8 Data")%, "1/16 Data")%, "1/64 Data")

%% Feed forward
% Preparing data and getting resonance
M_estimate=0.03127;

abs_H = abs(H);
[temp, index] = max(abs_H(30:100));
resonance_frequency = hz(30+index)

% Doing fit for mass (and spring constant)
ft = fittype('abs(1/(M*s^2))','dependent',{'H'},'independent',{'s'},'coefficients',{'M'});
f = fit(hz(2:30), abs_H(2:30), ft,'StartPoint', [M_estimate])
M=coeffvalues(f);
k_estimate=M*resonance_frequency^2

abs_H_rb = hz.^-2/(M);
Q = abs(H(index+30) - abs_H_rb(index+30));
d_estimate = M/(Q)%M*resonance_frequency/Q
hold on

% Do fit for damping
ft = fittype('abs(1./(M.*(i.*s)^2)-1./(M.*(i.*s)^2+d.*(i.*s)+k))','dependent',{'H'},'independent',{'s'},'coefficients',{'M','d','k'});
f = fit(hz(5:end), abs_H(5:end), ft,'StartPoint', [M,d_estimate, k_estimate])
coef = coeffvalues(f);
M=coef(1);
d=coef(2);
k=coef(3);

% Making a bode plot for verification
H_fit = 1./(M.*(hz.*1i).^2) - 1./(M.*(hz.*1i).^2+d.*(hz.*1i)+k);

H_mag_fit = abs(H_fit);
H_db_fit = mag2db(H_mag_fit);
H_angle_rad_fit = angle(H_fit);
H_angle_deg_fit = rad2deg(H_angle_rad_fit);

subplot(2,1,1)
semilogx(hz,H_db)
hold on
semilogx(hz,H_db_fit)%, "x")
title("Bode plot of FRF")
ylabel("Gain (dB)")
legend("FRF", "Fit")

subplot(2,1,2) 
semilogx(hz,H_angle_deg)
hold on
semilogx(hz,H_angle_deg_fit)%, "x")
ylabel("Phase (Degrees)")
xlabel("Angular frequency (s^{-1})")
legend("FRF", "Fit")
hold off

%% PSD analysis
nfft=250%250/4000%resolutie = fs/nfft
noverlap = nfft/2;
window = rectwin(length(e));
 
%subplot(2,1,1)
[apsd, hz_e] = periodogram(e, window, nfft, fs);%periodogram(e, window, nfft, fs)%pwelch(e, window, noverlap, nfft, 500);%periodogram(e_mask, window)
semilogx(hz_e, apsd)
xlabel("Frequency (Hz)")
ylabel("PSD")

r_dif=diff(r);
mask = abs(r_dif)>mean(abs(r_dif));
y_mask=y(mask);
%u_mask=u(mask);
e_mask=e(mask);

window_masked = rectwin(length(e_mask));%hann(nfft) of hamming()
figure
%subplot(2,1,2)
[apsd, hz_e] = periodogram(e_mask, window_masked, nfft, fs);%periodogram(e, window, nfft, fs)%pwelch(e, window, noverlap, nfft, 500);%periodogram(e_mask, window)
semilogx(hz_e, apsd)
xlabel("Frequency (Hz)")
ylabel("PSD with mask")