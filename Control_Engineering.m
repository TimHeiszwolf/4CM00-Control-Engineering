%% Measuring
%%% Settings
%p = raspi('10.55.0.1','pi','sperte123')
samples_start = 12*4e3;%4e3;
N_samples = 120*4e3;%4e3;% How many samples to measure. T = N_sampels/fs -> T*fs = N_sampels
noise_variance = 0.3;% Variance of the noise
N_signals = 2;% How many signals doclear you need to record?
ModelName = 'Reference_model_Control_Engineering';% Model should be Reference_model_Control_Engineering or change it if you renamed the model

N_samples = N_samples + samples_start;

%%% Runing measurment
% Check if the model is running, makes no sense to measure otherwise.
if ~strcmp(get_param(ModelName,'SimulationStatus'),'external')
   warning('The model is not running, can not perform a measurement.');
   return
end
% Get the noise block
BlockPaths = find_system(ModelName);
% Get the correct string
Index_noise = find(contains(BlockPaths,'Noise'));
% If we have more than 1 index, this means duplicate blocks, should not
% happen! Give a warning, and then try to trigger the first.
if length(Index_noise) > 1
    warning('There appear to be duplicate noise blocks in your model. Please remove all duplicate blocks.');
    Index_noise = Index_noise(1);
end

fprintf(['Starting measurement with white noise on the disturbance.\n',...
    'If you are planning on doing an FRF measurement, \n',...
    'think carefully about the signals you need and consider open/closed loop.\n']);

% Set the variance of the noise
set_param(BlockPaths{Index_noise},'variance', string(noise_variance));
pause(1.1);

% Do the measurement
measurement = SPERTE_Measure_And_Collect(p, N_signals, N_samples, ModelName);
pause(0.1);

% Turn of the noise
set_param(BlockPaths{Index_noise},'variance','0')

% Getting data
%r = measurement(samples_start+1:end,1);
u = measurement(samples_start+1:end,1);
%d = measurement(samples_start+1:end:,1);
y = measurement(samples_start+1:end,2);
%r = measurement(samples_start+1:end:,3);

N_samples = N_samples - samples_start;

%mdlWks = get_param(ModelName,'ModelWorkspace');%%%%% TODO TEST DIT
fs = 4000;%getVariable(mdlWks,'Fs');

%%% Plot the measurement
time = 1/fs * linspace(0, N_samples-1, N_samples);
figure(1); clf(1);
stairs(time, measurement(samples_start+1:end,1)); hold on; stairs(time, measurement(samples_start+1:end,2));
xlabel('Time [s]'); ylabel('Signal [volt,rad]');

% Saving data
[year, month, day] = ymd(datetime);
[hour, minute, second] = hms(datetime);
save("measurements/Measurement "+string(year)+"-"+string(month)+"-"+string(day)+"_"+string(hour)+"."+string(minute)+"."+string(round(second))+".mat")

%% Measurement from scope
samples_start = 12*0.25e3;

time = ScopeData(samples_start+1:end, 1);
dt = time(2)-time(1);
fs = 1/dt;
r = ScopeData(samples_start+1:end, 2);
y = ScopeData(samples_start+1:end, 3);
e = ScopeData(samples_start+1:end, 4);

% Saving data
[year, month, day] = ymd(datetime);
[hour, minute, second] = hms(datetime);
save("measurements/Measurement "+string(year)+"-"+string(month)+"-"+string(day)+"_"+string(hour)+"."+string(minute)+"."+string(round(second))+".mat")

%% Error calculations
all_data=false;% if false it only analyses the (first) period, else it does the entire dataset.
indexPeriod=2*12*fs;
y_axis_limit=0.2;
avg_length=100;

r_dif=diff(r);

% Removing turn around
mask = abs(r_dif)>mean(abs(r_dif));
mask_indexPeriod= abs(r_dif(1:indexPeriod))>mean(abs(r_dif(1:indexPeriod)));

% Checking error
%e=y-r;% Commented out due to e already being recorded properly.

adjusted_r=(0.5*r)*(0.8*max(e)/max(r));% A reference adjust to scale to the error. Can be used as overlay to give idea of reference signal.
adjusted_y=(0.5*y)*(0.8*max(e)/max(y));

turnaround_top=1000*(62.5-[max(y(1:indexPeriod)), max(y)])% Should not be more than 120+5
turnaround_bottom=1000*(62.5+[min(y(1:indexPeriod)), min(y)])% Should not be less than 0

max_error = 1000*[max(abs(e(1:indexPeriod))), max(abs(e))]
max_error_mask = 1000*[max(abs(e(mask_indexPeriod))), max(abs(e(mask)))]
rms_error = 1000*[rms(e(1:indexPeriod)), rms(e)]
rms_error_mask = 1000*[rms(e(mask_indexPeriod)), rms(e(mask))]

rms_error_mask_moving = zeros(1,length(e(mask)));
for i=avg_length:(length(e(mask))-avg_length)
    rms_error_mask_moving(i)=rms(e(i:(i+avg_length)));
end

% Making plots
if all_data
    plot(time, e)
    hold on
    plot(time, movmean(e, 0.25*fs), 'LineWidth', 2)
    plot(time, adjusted_r, 'LineWidth', 2)
    plot(time(avg_length:(length(e(mask))+avg_length-1)), rms_error_mask_moving, 'LineWidth', 2)
    axis([min(time) max(time) -y_axis_limit y_axis_limit])
else
    plot(time(1:indexPeriod), e(1:indexPeriod))
    hold on
    plot(time(1:indexPeriod), movmean(e(1:indexPeriod), 0.25*fs), 'LineWidth', 2)
    plot(time(1:indexPeriod), adjusted_r(1:indexPeriod), 'LineWidth', 2)
    plot(time(2*avg_length:(indexPeriod)), rms_error_mask_moving(avg_length:indexPeriod-avg_length), 'LineWidth', 2)
    axis([min(time(1:indexPeriod)) max(time(1:indexPeriod)) -y_axis_limit y_axis_limit])
end
%title("Error plot")
ylabel("Error (rad)")
xlabel("Time (s)")
legend("Error", "Moving mean of error", "Shape of reference", "Moving RMSE")
hold off
shg% Imediatly show the plot instead of waiting on performance region check.

%% Checking performance region
disp("Now checking constant velocity of reference.")
for i=1:length(r_dif)
    if (r(i)>-60 & r(i)<60) & abs(r_dif(i))<mean(abs(r_dif))%(r(i)>2.5 & r(i)<120)
        disp("Possible issue with constant velocity at r="+string(r(i))+" at index="+string(i))
    end
    if not(all_data) & i>indexPeriod
        disp("Only checked part of measurement")
        break
    end
end
disp("Finished constant velocity check.")

r_ft=fft(r);
f = (0:length(r)-1)*fs/length(r);
[temp, i]=max(r_ft(5:200));
period_reference=1.0/f(i+4)% Should be approximatly equal to the tracking period
%plot(f(1:100),abs(fft_r(1:100)))

