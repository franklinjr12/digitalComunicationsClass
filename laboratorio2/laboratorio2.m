close all;

%============================================================
%Geração do sinal y
%============================================================

Fsampling = 50;
Tsymbol = 1;
Vamp = 1;
Nbits = 10;

t = 0:1/Fsampling:Nbits-1/Fsampling;
bit1 = Vamp*ones(1,Fsampling);
bit0 = -Vamp*ones(1,Fsampling);
y = zeros(1,length(t));
for i=1:length(t)/Nbits:length(t)
    if (randi(100,1) <= 50)
        y(i:i+length(bit1)-1) = bit0;
    else
        y(i:i+length(bit1)-1) = bit1;
    end
end

figure;
plot(t,y);
xlabel('tempo (s)');
title('generator');
axis([0 Nbits -1.2*Vamp 1.2*Vamp]);

%============================================================
%1 - Transmissão Antipodal
%============================================================
h=[ones(1,Fsampling)];
r=conv(y,h)/Fsampling;
fim=length(r);
t=0:1/Fsampling:fim/Fsampling-1/Fsampling;
figure;
plot(t,r);
t_amostra=[Fsampling:Fsampling:Nbits*Fsampling-Fsampling];
r_amostra=r(t_amostra);

t_amostra=t_amostra/Fsampling-1/Fsampling;
hold
stem(t_amostra,r_amostra,'r')
xlabel('tempo (s)');
legend('Saída do Filtro','Saída Amostrada');
axis([0 max(t) -1.2*Vamp 1.2*Vamp])
%============================================================
%O meu critério do bit seria exatamente o valor de r_amostra
%============================================================

%============================================================
%2 - Adição do Ruído
%============================================================
noise_power = 1;
noise = wgn(1,length(y), noise_power);
figure;
plot(noise);
title('gaussian noise');

h=[ones(1,Fsampling)];
r=conv(y+noise,h)/Fsampling;
fim=length(r);
t=0:1/Fsampling:fim/Fsampling-1/Fsampling;
figure;
plot(t,r);
t_amostra=[Fsampling:Fsampling:Nbits*Fsampling-Fsampling];
r_amostra=r(t_amostra);

t_amostra=t_amostra/Fsampling-1/Fsampling;
hold
stem(t_amostra,r_amostra,'r')
xlabel('tempo (s)');
legend('Saída do Filtro','Saída Amostrada');
axis([0 max(t) -1.5*Vamp 1.5*Vamp])

%============================================================
%3 - Generalização
%============================================================
bit_error_rate = [];
wave_energy = [];
k=4;
max_bit_err = 10^-k; 
Nbits=100*10^k;
% Nbits=20;

Fsampling = 50;
Tsymbol = 1;
Vamp = 0.1;

y = zeros(1,length(t));

%iterate until 10dB
final_wave_energy_db = 10;
actual_wave_energy_db = 10*log10(Vamp^2*Fsampling*Tsymbol);
while(10*log10((Vamp^2)*Fsampling*Tsymbol) < final_wave_energy_db)

    %generate y signal with random bits
    t = 0:1/Fsampling:Nbits-1/Fsampling;
    bit1 = Vamp*ones(1,Fsampling);
    bit0 = -Vamp*ones(1,Fsampling);
    for i=1:length(t)/Nbits:length(t)
        if (randi(100,1) <= 50)
            y(ceil(i:i+length(bit1)-1)) = bit0;
        else
            y(ceil(i:i+length(bit1)-1)) = bit1;
        end
    end 

%     figure;
%     plot(t,y);
%     title('generated');
%     axis([0 max(t) -1.5*Vamp 1.5*Vamp]);
    
    %add transmission noise
    noise_power = 1;
    noise = wgn(1,length(y), noise_power);
    yr = y+noise;

    h=[ones(1,Fsampling)];
    r=conv(y+yr,h)/Fsampling;  
    fim=length(r);
    t=0:1/Fsampling:fim/Fsampling-1/Fsampling;
    
%     figure;
%     plot(t,r);
%     hold on;  
    
    %filter received signal
    h=[ones(1,Fsampling)];
    r=conv(yr,h)/Fsampling;
    fim=length(r);
    t=0:1/Fsampling:fim/Fsampling-1/Fsampling;
    t_amostra=[Fsampling:Fsampling:Nbits*Fsampling];
    r_amostra=r(t_amostra);

%     stem(t_amostra/Fsampling, r_amostra);
%     disp('pause');
%     pause;
    
    %calculate number of wrong bits
    n_wrong_bits = 0;
    aux = 1;
    step_value = length(y)/Nbits;
    for i = 1:step_value:length(y)
        generated_bit = mean(y(ceil(i:i+step_value-1)));
        if generated_bit > 0 
            generated_bit = Vamp;
        else
            generated_bit = -Vamp;
        end
        received_bit = r_amostra(aux);
        aux = aux+1;
        if received_bit > 0 
            received_bit = Vamp;
        else
            received_bit = -Vamp;
        end
%         generated_bit ~= received_bit
        if generated_bit ~= received_bit
            n_wrong_bits = n_wrong_bits+1;
        end
    end

%     disp('pause');
%     pause;
%     close all;
    
    %calculate statistics for making final graphs
    actual_wave_energy_db = 10*log10(Vamp^2*Fsampling*Tsymbol)
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

%============================================================
%O resultado acabou sendo bem parecido com o do laboratório anterior
%porém vale ressaltar que no passado foi usado um filtro passa baixas
%e o ruído foi gerado com uma função rand e não o gaussiano
%============================================================

disp('pause');
pause;

%============================================================
%4 - 4-PAM
%============================================================

bit_error_rate = [];
wave_energy = [];
k=4;
max_bit_err = 10^-k; 
Nbits=100*10^k;
% Nbits=20;
M=4;

Fsampling = 50;
Tsymbol = 1;
Vamp = 0.1;

y = zeros(1,length(t));

%iterate until 10dB
final_wave_energy_db = 10;
actual_wave_energy_db = 10*log10(Vamp^2*Fsampling*Tsymbol);
while(10*log10((Vamp^2)*Fsampling*Tsymbol) < final_wave_energy_db)

    %generate y signal with random bits
    t = 0:1/Fsampling:Nbits-1/Fsampling;
    bit3 = Vamp*ones(1,Fsampling);
    bit2 = -Vamp*ones(1,Fsampling);
    bit1 = Vamp/2*ones(1,Fsampling);
    bit0 = -Vamp/2*ones(1,Fsampling);    
    for i=1:length(t)/Nbits:length(t)
        if (randi(100,1) <= 25)
            y(ceil(i:i+length(bit1)-1)) = bit0;
        elseif (randi(100,1) <= 50)
            y(ceil(i:i+length(bit1)-1)) = bit1;
        elseif (randi(100,1) <= 75)
            y(ceil(i:i+length(bit1)-1)) = bit2;
        else
            y(ceil(i:i+length(bit1)-1)) = bit3;
        end
    end 

%     figure;
%     plot(t,y);
%     title('generated');
%     axis([0 max(t) -1.5*Vamp 1.5*Vamp]);
    
    %add transmission noise
    noise_power = 1;
    noise = wgn(1,length(y), noise_power);
    yr = y+noise;

    h=[ones(1,Fsampling)];
    r=conv(y+yr,h)/Fsampling;  
    fim=length(r);
    t=0:1/Fsampling:fim/Fsampling-1/Fsampling;
    
%     figure;
%     plot(t,r);
%     hold on;  
    
    %filter received signal
    h=[ones(1,Fsampling)];
    r=conv(yr,h)/Fsampling;
    fim=length(r);
    t=0:1/Fsampling:fim/Fsampling-1/Fsampling;
    t_amostra=[Fsampling:Fsampling:Nbits*Fsampling];
    r_amostra=r(t_amostra);

%     stem(t_amostra/Fsampling, r_amostra);
%     disp('pause');
%     pause;
    
    %calculate number of wrong bits
    n_wrong_bits = 0;
    aux = 1;
    step_value = length(y)/Nbits;
    for i = 1:step_value:length(y)
        generated_bit = mean(y(ceil(i:i+step_value-1)));
        if generated_bit < -Vamp/2 
            generated_bit = -Vamp;
        elseif generated_bit < 0
            generated_bit = -Vamp/2;
        elseif generated_bit < Vamp
            generated_bit = Vamp/2;
        else
            generated_bit = Vamp;
        end
        received_bit = r_amostra(aux);
        aux = aux+1;
        if received_bit < -Vamp/2 
            received_bit = -Vamp;
        elseif received_bit < 0
            received_bit = -Vamp/2;
        elseif received_bit < Vamp
            received_bit = Vamp/2;
        else
            received_bit = Vamp;
        end
%         generated_bit ~= received_bit
        if generated_bit ~= received_bit
            n_wrong_bits = n_wrong_bits+1;
        end
    end

%     disp('pause');
%     pause;
%     close all;
    
    %calculate statistics for making final graphs
    actual_wave_energy_db = 10*log10(Vamp^2*Fsampling*Tsymbol)
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

%============================================================
%Mesmo sem ter alterado a amplitude da onda o resultado não piorou
%tanto quanto eu pensei que iria piorar. O erro aumentou 3 vezes
%============================================================
