%Programa para calcular un filtro pasa bandas
%Aproximacion butterworth y una realizacion Silent key
clear all;
close all;
clc;

%Paso 1 Especificaciones del filtro

Ap = 1;
As1 = 35;
As2 = 40;
Fs1 = 1000;
Fc1 =2000;
Fs2=8000;
Fc2= 4000;

%Paso 2 Normalizacion BPF
Apn = Ap;
Asn = max([As1, As2]);

if((Fc1/Fs1)>(Fs2/Fc2))
  Fs1_a = (Fc1 * Fc2)/Fs2;
  Fs2_a = Fs2;
else
  Fs1_a = Fs1;
  Fs2_a = (Fc1 * Fc2)/Fs1;
end

Fsn = (Fs2_a - Fs1_a)/(Fc2-Fc1);


% PASO 3. Calculo del orden del filtro pasa bajas normalizado
epsilon=sqrt(10^(Apn/10)-1);
N=log10((10^(Asn/10)-1)/epsilon^2)/(2*log10(Fsn));
N=ceil(N);

% PASO 4. Calculo de las bicuadraticas
if (mod(N,2))
% N es impar
  Num_N(1,:)=[0 0 epsilon^(-1/N)];
  Den_N(1,:)=[0 1 epsilon^(-1/N)];
  for k=1:(N-1)/2
    Num_N(k+1,:)=[0 0 (epsilon^(-1/N)^2)];
    Den_N(k+1,:)=[1 -2*epsilon^(-1/N)*cos(pi*(N+2*k-1)/(2*N)) (epsilon^(-1/N)^2)];
  endfor
else
% N es par
  for k=1:N/2
    Num_N(k,:)=[0 0 (epsilon^(-1/N)^2)];
    Den_N(k,:)=[1 -2*epsilon^(-1/N)*cos(pi*(N+2*k-1)/(2*N)) (epsilon^(-1/N)^2)];
  endfor
end

% PASO 5. Desnormalizacion del filtro
a = 4*Fc1*Fc2*pi^2;
b = 2*pi*(Fc2-Fc1);
if (mod(N,2))
##% N es impar
  Num(1,:)=[0 Num_N(1,3)*b 0];
  Den(1,:)=[1 Den_N(1,3)*b a];
  for k=1:(N-1)/2
    p = [1 Den_N(k+1,2)*b (Den_N(k+1,3)*b^2 + 2*a) Den_N(k+1,2)*a*b a^2];
    raices = roots(p);
    Num(2*k,:) = [0 sqrt(Num_N(k+1,3))*b 0];
    Den(2*k,:) = [1 -2*real(raices(1)) abs(raices(1))^2]
    Num(2*k+1,:) = [0 sqrt(Num_N(k+1,3))*b 0];
    Den(2*k+1,:) = [1 -2*real(raices(3)) abs(raices(3))^2]
  endfor
else
% N es par
  for k=1:N/2
    p = [1 Den_N(k,2)*b (Den_N(k,3)*b^2 + 2*a) Den_N(k,2)*a*b a^2];
    raices = roots(p);
    Num(2*k-1,:) = [0 sqrt(Num_N(k,3))*b 0];
    Den(2*k-1,:) = [1 -2*real(raices(1)) abs(raices(1))^2]
    Num(2*k,:) = [0 sqrt(Num_N(k,3))*b 0];
    Den(2*k,:) = [1 -2*real(raices(3)) abs(raices(3))^2]
  endfor
end


% PASO 6. Verificacion de la respuesta en frecuencia
f=logspace(1,5,1000);
s=1j*2*pi*f;

H=1;
  for k=1:N
    H=H.*((Num(k,1)*s.^2 + Num(k,2).*s + Num(k,3))./(Den(k,1)*s.^2 + Den(k,2).*s + Den(k,3)));
  endfor
figure;
subplot(2,1,1); semilogx(f,20*log10(abs(H))); grid on; hold on;
                semilogx([10 Fs1 Fs1],[-As1 -As1 0],'r');
                semilogx([Fc1 Fc1 Fc2 Fc2],[-100 -Ap -Ap -100],'r');
                semilogx([Fs2 Fs2 10^5],[0 -As2 -As2], 'r');
                axis([10 10^5 -100 10]);
subplot(2,1,2); semilogx(f,angle(H)); grid on;

