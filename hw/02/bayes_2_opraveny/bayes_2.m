%% Bayesiánska analýza - úloha číslo 2
% spracovali Petr Kukuczka, Monika Dvořáčková,
% Patrik Suchánek, Michaela Kozaňáková

clear all
close all
clc

addpath('.\Support')
load('data.mat')

U_1 = DATA(:,1);
HICP_1 = DATA(:,2);

U_0 = [U_1(1:71)];
U_1 = [U_1(2:72)];

U_diff = [U_1-U_0];
HICP_0 = [HICP_1(1:71)];
HICP_1 = [HICP_1(2:72)];

%% 1. Definovanie dát a apriorných hyperparametrov
 b0 =[1.3; 0.7; -0.25; -0.5];       
%zadanie apriornej strednej hodnoty pre vektor parametra b
 V0 = diag([1^2;0.5^2;0.5^2;1^2]); 
%zadanie apriornej kovariančnej matice parametra b
 v0 = 20;                          
%zadanie apriorneho počtu stupňov voľnosti
 s2_0 = 10^2;                       
%apriorni rozptyl náhodných zložiek
 h0 = 1/s2_0;                       
%vygenerovanie apriornej presnosti chyby
y = HICP_1;  %vysvetľovaná premenná
X = [ones(size(y)) HICP_0 U_1 U_diff]; %vysvetľujúca premenná
[N,k] = size(X);

%% 2. Inicializácia Gibbsova vzorkovača
 S1 = 50000;    %počet ponechaných vzorkov
 S0 = 50000+1;  %počet vyhodených vzorkov
 S = S1+S0+1;   %celkový počet generovaných vzorkov+počiatočná hodnota
 
 
 b = zeros(k,S);%priestor pre ukladanie vzorku b
 h = zeros(1,S);%vzorky pro h
 
%časť pre nastavenie počatočných hodnôt 
b(:,1) = [0;0;0;0];
h(:,1) = h0;
 
%% Gibbsuv vzorkovač
eta = 1.5; %táto hodnota je zadaná ľubovoľne

if ((eta >= 1) || (eta <= 0))
   for s=2:S
    %a) podmienená hustota p(b|h,y)
    %b~N(b1,V1) normálne rozdelenie
    V1 = inv(inv(V0)+h(1,s-1)*(X'*X));
    b1 = V1*(inv(V0)*b0+h(1,s-1)*(X'*y)); 
    b(:,s) = norm_rnd(V1)+b1;
        
    %b) podmienená hustota p(h|b,y)
    %h~G(h1,v1)rozdelenie pre h
    v1 = N+v0;
    h1 = (1/v1*((y-X*b(:,s))'*(y-X*b(:,s))+v0*h0^-1))^-1;
    h(s) = gamm_rnd_Koop(h1,v1,1);      
   end   
    
% Posteriorné momenty
% Apriorné a posteriorné stredné hodnoty
% Vyhodnotenie prvých S0+1 pozorovaní
  b(:,1:S0+1)=[];
  h(1:S0+1) =[];
 
  E_b = mean(b,2);
  E_h = mean(h);
  D_b = var(b,1,2);
  D_h = var(h,1); 
 
% Prevedenie výpočtu štrukturálnych parametrov
  alpha = E_b(2);
  beta = E_b(3) + E_b(4);
  Z = E_b(1)/(-beta);  
  eta = E_b(4)/beta;
end

U = (eta*U_0) + Z;
 
%% Grafické prevedenie posteriorných hustôt
 figure
  for ii=1:k
   subplot(2,2,ii)
   hist(b(ii,:),20)
   xlabel(['\lambda_',num2str(ii)])
  end
 
 saveas(gcf,'Obrazek1.png');
  
 figure
  hist(h,20)
  xlabel('h')
  saveas(gcf,'Obrazek2.png');
 
%Grafický pohľad na konvergenciu lambda
 figure
  for ii=1:k
   subplot(2,2,ii)
   plot(b(ii,1:500:end))
   xlabel(['\lambda_',num2str(ii)])
  end
 saveas(gcf,'Obrazek3.png');
  
%% Konvergenčná diagnostika
 CD_b = Geweke(b');
 CD_h = Geweke(h');

%% %Savage-Dickey pomer hustôt obecne pre obmedzenie R*lambda = r 
SD = zeros(4,1);
SD_nom = zeros(S1,1);

%zapísanie obmedzenia
R = [1 0 0 0];
r = 0;

for s=1:S1
    V1 = inv(inv(V0)+h(1,s)*(X'*X));
    b1 = V1*(inv(V0)*b0+h(1,s)*(X'*y)); 
    %podmienená hustota pre b|h~N(b1,V1) => R*b|h~N(R*b1;R*V1*R')
    SD_nom(s) = 1/(2*pi)^(1/2)*(det(R*V1*R'))^(-1/2)*exp(-1/2*(r-R*b1)'*(R*V1*R')^-1*(r-R*b1));   
end

%Časť pre vyhodnotenie menovateľa
%podmienená hustota pre b|h~N(b0,V0) => R*b|h~N(R*b0;R*V0*R')
SD_denom = 1/(2*pi)^(1/2)*(det(R*V0*R'))^(-1/2)*exp(-1/2*(r-R*b0)'*(R*V0*R')^-1*(r-R*b0));   
%S-D pomer hustôt
SD(1,1) = mean(SD_nom)/SD_denom;

%Zapísanie obmedzenia
R = [0 1 0 0];
r = 0;

for s=1:S1
    V1 = inv(inv(V0)+h(1,s)*(X'*X)); 
    b1 = V1*(inv(V0)*b0+h(1,s)*(X'*y)); 
    %podmienená hustota pre b|h~N(b1,V1) => R*b|h~N(R*b1;R*V1*R')
    SD_nom(s) = 1/(2*pi)^(1/2)*(det(R*V1*R'))^(-1/2)*exp(-1/2*(r-R*b1)'*(R*V1*R')^-1*(r-R*b1));   
end

%vyhodnotenie menovateľa
%podmienená hustota pre b|h~N(b0,V0) => R*b|h~N(R*b0;R*V0*R')
SD_denom = 1/(2*pi)^(1/2)*(det(R*V0*R'))^(-1/2)*exp(-1/2*(r-R*b0)'*(R*V0*R')^-1*(r-R*b0));   
%S-D pomer hustôt
SD(2,1) = mean(SD_nom)/SD_denom;

%Zapísanie obmedzenia 
R = [0 0 1 0];
r = 0;

for s=1:S1
    V1 = inv(inv(V0)+h(1,s)*(X'*X)); 
    b1 = V1*(inv(V0)*b0+h(1,s)*(X'*y)); 
    %podmienená hustota pre b|h~N(b1,V1) => R*b|h~N(R*b1;R*V1*R')
    SD_nom(s) = 1/(2*pi)^(1/2)*(det(R*V1*R'))^(-1/2)*exp(-1/2*(r-R*b1)'*(R*V1*R')^-1*(r-R*b1));   
end

%vyhodnotenie menovateľa
%podmienena hustota pre b|h~N(b0,V0) => R*b|h~N(R*b0;R*V0*R')
SD_denom = 1/(2*pi)^(1/2)*(det(R*V0*R'))^(-1/2)*exp(-1/2*(r-R*b0)'*(R*V0*R')^-1*(r-R*b0));   
%S-D pomer hustôt
SD(3,1) = mean(SD_nom)/SD_denom;

%Zapísanie obmedzenia 
R = [0 0 0 1];
r = 0;

for s=1:S1
    V1 = inv(inv(V0)+h(1,s)*(X'*X));
    b1 = V1*(inv(V0)*b0+h(1,s)*(X'*y));  
    %podmienená hustota pre b|h~N(b1,V1) => R*b|h~N(R*b1;R*V1*R')
    SD_nom(s) = 1/(2*pi)^(1/2)*(det(R*V1*R'))^(-1/2)*exp(-1/2*(r-R*b1)'*(R*V1*R')^-1*(r-R*b1));   
end
%vyhodnoceni jmenovatele
%podminena hustota b|h~N(b0,V0) => R*b|h~N(R*b0;R*V0*R')
SD_denom = 1/(2*pi)^(1/2)*(det(R*V0*R'))^(-1/2)*exp(-1/2*(r-R*b0)'*(R*V0*R')^-1*(r-R*b0)); 

%S-D pomer hustôt
SD(4,1) = mean(SD_nom)/SD_denom;

%% Časť pre výpočet posteriornich pravdepodobností modelu
p3 = 1/(1+SD(3,1)+SD(4,1));
p2 = SD(4,1)*p3;
p1 = SD(3,1)*p3;

%% Časť pre výpis výsledkov
fprintf('          Odhad parametrov lambda a h (v zátvorkách sú uvedené smerodatne odchylky)          \n');
fprintf('                  E Prior    E Posterior   Gewekeho CD     Post. Odds\n');
for i=1:k
fprintf('Lambda %2u          %2.4f       %2.4f      %3.4f         %3.4f \n',[i b0(i) E_b(i) CD_b.CD(i) SD(i)]);
fprintf('                   (%2.4f)     (%2.4f)       \n',[sqrt(V0(i,i)) sqrt(D_b(i))]);
end
fprintf('h                  %2.4f       %2.4f       %3.4f    \n',[h0 E_h CD_h.CD]);
fprintf('                   (%2.4f)    (%2.4f)       \n',[sqrt(2*h0^2/v0) sqrt(D_h)]);
fprintf('--------------------------------------------------------------------------------\n');
fprintf('          Štrukturálne parametry         \n');
fprintf('                    alpha     beta     eta         Z\n');
fprintf('Lucembursko         %2.4f    %2.4f   %2.4f   %2.4f \n',[alpha beta eta Z]);

fprintf('--------------------------------------------------------------------------------\n');
fprintf('          Posteriorné pravdepodobnosti modelu         \n');
fprintf('                     eta = 0     eta = 1     eta patri (0,1)     \n');
fprintf('Lucembursko          %2.4f      %2.4f         %2.4f   \n',[p2 p1 p3]);
 
%% Monte Carlo integrácia
%preddefinovanie matice NAIRU

NAIRU = NaN(length(U_1),S1);

%simulácia 
for ss = 1:S1
alpha = b(2,ss);
beta = b(3,ss) + b(4,ss);
Z = b(1,ss)/(-beta);
eta = b(4,ss)/beta;

NAIRU(:,ss) = eta*U_0+Z;
end


NAIRU_mean = median(NAIRU,2);
% HPDI
NAIRU_05 = quantile(NAIRU,0.05,2);
NAIRU_95 = quantile(NAIRU,0.95,2);

%% Vykreslenie trajektorie NAIRU
startDate = datenum('2000-01-01');
endDate = datenum('2017-12-31');

%Vytvorenie dátových bodov 'xdata' (osa x) zodpovedajúcich počtov mesiacov:
xData = linspace(endDate, startDate, 71);
figure
plot(xData, NAIRU_mean,'k', xData, U_1,'--', xData, NAIRU_05,':k', xData, NAIRU_95,':k')
 datetick('x','QQ yyyy','keepticks')

 %nastavenie limity osy x
xlim([startDate-30 endDate+30]); 

 title('Obrazek 1: NAIRU a nezamestnanost')
 legend('NAIRU','Nezamestnanost','95% HPDI','95% HPDI')
 saveas(gcf,'Obrazek4.png');
