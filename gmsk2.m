clear all
gmsk_mod = comm.GMSKModulator('BitInput', true);
gmsk_demod = comm.GMSKDemodulator('BitOutput', true);
awgn_channel = comm.AWGNChannel('NoiseMethod', ...
                    'Signal to noise ratio (SNR)','SNR',0);
error_rate = comm.ErrorRate('ReceiveDelay', gmsk_demod.TracebackDepth);
signal_size = 1024;
iterations = 1000;
for i = 1:iterations
    my_signal = randi([0,1],signal_size,1)';
    my_signal_b = my_signal;
    signal_gmsk_mod = step(gmsk_mod, my_signal_b(:));
    signal_gmsk_noisy = step(awgn_channel, signal_gmsk_mod);
    signal_gmsk_demod = step(gmsk_demod, signal_gmsk_noisy);
    error_stat = step(error_rate, my_signal_b(:), signal_gmsk_demod);
end
fprintf('Error rate = %f\nNumber of errors = %d\n', ...
      error_stat(1), error_stat(2))