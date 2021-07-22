close all;

% Respostas:
% O sinal original é recuperado?
% Sim!
% A BER muda em relação ao programa original quando introduzimos ruído? Porque?
% Sim! A BER foi substancialmente menor em relação ao programa original.
% E como fica a questão do atraso do filtro para os diferentes fatores de excesso de faixa?
% Quanto menor o fator mais atraso o filtro tem.
% Finalmente, porque não usamos na prática o menor valor de excesso de faixa possível para o raiz de cosseno levantado?
% Porque isso torna inviável o uso do filtro devido um atraso grande na
% reposta.

%initial values and arrays

Fs = 24000;
T = 1/8000;
r = Fs*T; % Fator de oversampling
t = -5*T:1/Fs:5*T;
t = t + 1e-10;
alfa = 0.5;
% alfa = 0.00001;

bit_error_rate = [];
wave_energy = [];
k=4;
max_bit_err = 10^-k;
Nbits=100*10^k;

%generate sqrt of upper cosine
h = rcosine(1/T,Fs,'sqrt',alfa,5);

%generate bit sequence
bits = zeros(1,Nbits);
for i=1:Nbits
    bits(i) = rand > 0.5;
end
s=2*bits-1;
s_up = zeros(1,length(s)*r);

%cycle until db threshold
Vamp = 0.1;
final_wave_energy_db = 10;
actual_wave_energy_db = 10*log10(Vamp^2*Fs*T);
while(10*log10((Vamp^2)*Fs*T) < final_wave_energy_db)
    
    % sequência com oversampling
    s_up(1:r:r*length(s)) = s*Vamp;
    x=conv(s_up,h);
    
    %add transmission noise
    noise = (rand(1,length(x))-0.5)*2;
    s_trans = x+noise;
        
    %receive signal
    receiver = conv(s_trans,h);
    eixo = (0:length(receiver)-1)*1/Fs;
    
    bits_received = zeros(1,Nbits);
    j=1;
    t_start_symbol = 0.00125;
    for i=1:length(receiver)
        if (eixo(i) >= t_start_symbol && eixo(i) < T*Nbits + t_start_symbol ...
                && j <= Nbits)            
            if (mod(eixo(i),T)==0)
                if (receiver(i) > 0)
                    bits_received(j)=1;
                    j=j+1;
                else
                    bits_received(j)=0;
                    j=j+1;
                end
            end
        end
    end
           
    %calculate statistics for making final graphs
    n_wrong_bits = sum( bits ~= bits_received);
    actual_wave_energy_db = 10*log10(Vamp^2*Fs*T)
    bit_error_rate = [bit_error_rate, (n_wrong_bits/Nbits)];
    wave_energy = [wave_energy , actual_wave_energy_db];
    
    %update Amplitude for next iteraction
    Vamp = Vamp+0.1;
end

figure();
semilogy(wave_energy, bit_error_rate);
axis([0 max(wave_energy) 0 1.2*max(bit_error_rate)]);
xlabel('wave energy (dB)');
ylabel('bit error rate (log)');
disp('Final amplitude: ');
disp(Vamp);