%d -> xp es la señal de entrada
%y -> x_hat es la señal de salida del filtro
%x -> r la referencia 
% X -> R ruido refere
% c -> coeficientes
clc
close all;
clear all;  
 
%--------------------------------------------------------------------------
%                      Se cargan las señales de audio
%--------------------------------------------------------------------------

[senal,FS]= audioread('its-me-mario.mp3'); 
[ruido,FS]= audioread('marioSoloRuido.wav');
[senal_ruido,FS]= audioread('marioConRuido.wav');
x = ruido(1:size(ruido),1); % Señal de Referencia
d = senal_ruido; % Señal deseada 
    
%--------------------------------------------------------------------------
%                      Parámetros del Filtro
%--------------------------------------------------------------------------
% Calculamos el numero de muestras
NMuest=size(d);       % setea el nº de muestras de la señal
Ncoef= 1024;          %orden del filtro adaptativo, cantidad de coeficientes
W = zeros (Ncoef,1);  % genera una matriz cuyos elementos son todos ceros
Mu= 0.005;            % set the step-size constant
%inicializamos la primera colñumna del LMS
%potencia = sum(r.^2); % potencia del ruido
%calidad = 1; 
%Mu = ((2*calidad)/potencia)/3 
          
%--------------------------------------------------------------------------
%                         Algoritmo LMS 
%--------------------------------------------------------------------------
for n=Ncoef:NMuest 
	X = x(n:-1:n-Ncoef+1); 
	y(n) = W'*X;           % Salida del filtro, convolucion de W y X
	e(n) = d(n)-y(n);      % calculamos el error
	W = W+Mu* e(n)*X;      % actualizamos los coeficientes // 2u saber xq se elige Mu
    error(n) =(e(n)^2);
    MSE(n-Ncoef+1)= 10*log10(error(n));  
end;    
    
%--------------------------------------------------------------------------
%                         Guardar sonidos 
%--------------------------------------------------------------------------

audiowrite('C:\Users\rafag\Desktop\TP\filtrado.wav', e, FS);

sound(e,FS) % Escuchar salida del Filtro
%sound(d,FS) % Escuchar Señal+ruido
%sound(x,FS) % Escuchar Ruido
 
%--------------------------------------------------------------------------
%                         Calculos del FFT
%--------------------------------------------------------------------------
L=length(e); %Filtro
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
fft_signal=fft(e,NFFT)/L;
f = FS/2*linspace(0,1,NFFT/2);

L2=length(d);
NFFT_2 = 2^nextpow2(L2); % Next power of 2 from length of y
fft_signal2=fft(d,NFFT_2)/L2;
f2 = FS/2*linspace(0,1,NFFT_2/2); 

L3=length(senal);%limpia senal
NFFT_3 = 2^nextpow2(L3); % Next power of 2 from length of y
fft_signal3=fft(senal,NFFT_3)/L3;
f3 = FS/2*linspace(0,1,NFFT_3/2);

%--------------------------------------------------------------------------
%                         Prueba  lms del matlab
%--------------------------------------------------------------------------
if size(x, 2) == 2
    x = mean(x, 2); % Convierte a mono
end
if size(d, 2) == 2
    d = mean(d, 2); % Convierte a mono
end
ha = dsp.LMSFilter('Length', 512, 'StepSize', 0.001);
[y, e1] = step(ha, x, d);
legend('Salida Filtro calculado por LMS');
for n=1:NMuest
    MSE_matlab(n)= (e1(n)^2);
end  

   
%--------------------------------------------------------------------------
%                         Grafica del error
%--------------------------------------------------------------------------
plot((senal+x),'m');hold on;plot(d,'r');axis([0 100000 -2 2])
%subplot(2,2,3);plot(senal);subplot(2,2,4);plot(e);title('Rojo Ruido, Verde senal+ruido, Azul senal');
  
%--------------------------------------------------------------------------
%                         Grafica del MSE 
%--------------------------------------------------------------------------
figure(1)
semilogy([1:800],MSE(1:800),'g'); 
title('Curva de Convergencia');
xlabel('Nro de iteraciones');ylabel('MSE');grid on;  

%--------------------------------------------------------------------------
%                         Graficas del FFT
%--------------------------------------------------------------------------
figure(2)
plot(f2,2*abs(fft_signal2(1:NFFT_2/2)),'g') 
xlabel('Frequency (Hz)')
ylabel('|Y(f)|'); hold on; %axis([0 4000 0 1]);
% -------------
plot(f3,2*abs(fft_signal3(1:NFFT_3/2)),'m') %senal
xlabel('Frequency (Hz)')
ylabel('|Y(f)|');
% -------------
plot(f,2*abs(fft_signal(1:NFFT/2)),'r') 
title('ESPECTROS DE SEÑALES') 
xlabel('Frequency (Hz)')
ylabel('|Y(f)|');
axis([0 700 0 0.05]);legend('Senhal + Ruido','Senhal sin ruido','Senhal Filtro');grid on;

%--------------------------------------------------------------------------
%                  Grafica de comparacion  lms del matlab
%--------------------------------------------------------------------------
figure(3)
subplot(2,1,1)   
plot(e,'m') 
title('Curva de error');
hold on 
plot(e1); 
grid on; xlabel('Nro de iteraciones')
ylabel('Error'); 
legend('LMS diseñado','dsp.LMSFilter') % Ajustado para coincidir con la cantidad de series de datos

subplot(2,1,2)  
plot(MSE_matlab,'k');
hold on 
plot(MSE,'m'); % Asumiendo que MSE está correctamente calculado
xlabel('Nro de iteraciones'); 
title('Curva de Convergencia');
ylabel('MSE'); 
legend('LMS diseñado','dsp.LMSFilter'); % Ajustado para coincidir con la cantidad de series de datos
grid on;
