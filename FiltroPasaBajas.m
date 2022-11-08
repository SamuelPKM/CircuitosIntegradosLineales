% script para el calculo de un filtro pasa bajas mediante
% una aproximacion Butterworth y una realizacion Sallen-Key
clear all;
close all;
clc;

% PASO 1. Especificaciones del filtro:
Ap=1;
As=45;
Fc=1000;
Fs=4000;

R1_REF=2200;
R2_REF=2200;

% PASO 2. Normalizacion a un filtro pasa bajas
Apn=Ap;
Asn=As;
Fsn=Fs/Fc;

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
if (mod(N,2))
% N es impar
  Num(1,:)=[0 0 Num_N(1,3)*(2*pi*Fc)];
  Den(1,:)=[0 1 Den_N(1,3)*(2*pi*Fc)];
  for k=1:(N-1)/2
    Num(k+1,:)=[0 0 Num_N(k+1,3)*(2*pi*Fc)^2];
    Den(k+1,:)=[1 Den_N(k+1,2)*(2*pi*Fc) Den_N(k+1,3)*(2*pi*Fc)^2];
  endfor
else
% N es par
  for k=1:N/2
    Num(k,:)=[0 0 Num_N(k,3)*(2*pi*Fc)^2];
    Den(k,:)=[1 Den_N(k,2)*(2*pi*Fc) Den_N(k,3)*(2*pi*Fc)^2];
  endfor
end

% PASO 6. Verificacion de la respuesta en frecuencia
f=logspace(1,5,1000);
s=1j*2*pi*f;

H=1;
if (mod(N,2))
% N es impar
  H=H.*((Num(1,1)*s.^2 + Num(1,2).*s + Num(1,3))./(Den(1,1)*s.^2 + Den(1,2).*s + Den(1,3)));
  for k=1:(N-1)/2
    H=H.*((Num(k+1,1)*s.^2 + Num(k+1,2).*s + Num(k+1,3))./(Den(k+1,1)*s.^2 + Den(k+1,2).*s + Den(k+1,3)));
  endfor
else
% N es par
  for k=1:N/2
    H=H.*((Num(k,1)*s.^2 + Num(k,2).*s + Num(k,3))./(Den(k,1)*s.^2 + Den(k,2).*s + Den(k,3)));
  endfor
end
figure;
subplot(2,1,1); semilogx(f,20*log10(abs(H))); grid on; hold on;
                semilogx([10 Fc Fc],[-Ap -Ap -100],'r');
                semilogx([Fs Fs 10^5],[0 -As -As],'r');
                axis([10 10^5 -100 10]);
subplot(2,1,2); semilogx(f,angle(H)); grid on;

% PASO 6. Calculo de los compenetes de la estructura Sallen-Key
if (mod(N,2))
% El filtro es impar

else
% El filtro es par
  for k=1:N/2
    R(k,:)=[R1_REF R2_REF];
    C1=(1/R(k,1) + 1/R(k,2))/Den(k,2);
    C2=1/(Den(k,3)*R(k,1)*R(k,2)*C1);
    C(k,:)=[C1 C2];
  endfor
end





































