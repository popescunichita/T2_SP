P = 40;
D = 21;
N = 50;
W = 2*pi / P;   %Pulsatia unghiulara a semnalului dreptunghiular

pas = P / 100;  
t = (-2 * P) : pas : (2 * P);     %Pasul de esantionare

coefSFC = zeros(1,N);   %Initializarea vectorului de coeficienti SFC
coefSFA = zeros(1,N);   %Initializarea vectorului de coeficienti SFA

x = square( W * t, D); %Semnalul initial
x_init = @(t,k) square( W * t, D).*exp( -1j * k * W * t); %Semnal descompus in SFC
x_reconstr = 0;    %Initializarea semnalului reconstruit
coefxk = (1 / P) * integral(@(t) x_init(t,0),0,P);   %Coeficientii Xk ai semnalului


for k = 1:50    
    coefSFC(k) = (1 / P) * integral(@(t) x_init(t,k - 25 ),0,P);   
    x_reconstr = x_reconstr + coefSFC(k) * exp( 1j * (k - 25) * W * t);  %Reconstruirea semnalului
end                                                                               

figure(1);      %Reprezentarea semnalului dreptunghiular initial si reconstruit

plot(t, x_reconstr,'r', t, x,'b');
title('Semnalul dreptunghiular inital si reconstruit');
legend('semnal reconstruit','semnal initial');


coefSFA(1) = abs(coefxk);   %Coeficientii SFA - amplitudiniile din spectru
for k = 1:N
      coefSFA(k+1) = 2 * abs(coefSFC(k)); 
end

figure(2);      %Reprezentarea spectrului de amplitudini

stem([0:N] * W, coefSFA,'b'); 
title('Spectrul semnalului dreptunghiular');

%Conform teoriei seriilor Fourier, orice semnal periodic poate fi aproximat printr-o suma
%de sinusuri si cosinusuri de frecvente diferite si coeficienti diferiti, acesti coeficienti
%reprezentand spectrul de amplitudini al semnalului periodic initial. 
%Prin folosirea unui numar finit de termeni(N=50) semnalul reconstruit pe baza
%spectrului se apropie de cel initial avand insa o marja de eroare.