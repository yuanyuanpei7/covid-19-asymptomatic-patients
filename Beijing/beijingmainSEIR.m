,close all;clear all; clc;
 
num = xlsread('beijing.xlsx');
store=num(1:length(num(:,1)),[1,3,2]); %%%有症状感染者，治愈者，无症状感染者
Ipresent=[store(:,1),store(:,2),store(:,3)];%%%%%%from the 25 days start,the data is not correct at all before the 25th day;

%%% all togegher
%figure(4)
%plot(1:length(Ipresent(:,1)),Ipresent(:,1))

%%%%%%%%%%% investigate data
Ipresent=Ipresent(1:30,:); %%%second peak
%Ipresent=Ipresent(320:360,:); %%%third peak
%Ipresent=Ipresent(365:420,:); %%% fourth peak
%Ipresent(:,2)=Ipresent(:,2)-Ipresent(1,2)+Ipresent(1,1)/10; %%% cancle the first recover
t0=1:length(Ipresent(:,1));
 
N0=3e3; %total pupulation take into 作用群体人数需要修正

E0=Ipresent(1,1)+50;   %%%初值太低，需要修正
I0=Ipresent(1,1)+50;   %%%初值太低，需要修正
R0=Ipresent(1,2);
A0=Ipresent(1,3);
S0=N0-E0-I0-R0;

ini=[S0 E0 I0  A0 R0];
tspan=1:length(Ipresent(:,1));

seir=SEIR(t0,Ipresent,ini); %%% object
%%%%%%%%%%%%%%%%%%%% standard model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

c1=0.1;
c2=3;
Pi=0.5;
beta=0.25;
alpha1=0.2;
alpha2=0.2;
gama1=0.15;
gama2=0.15;
N=N0;
decays=0.1;

para=[c1 c2 Pi beta alpha1 alpha2 gama1 gama2 N0 decays];
lb=[0.01 0.9 0.2 0.1 0.16 0.16 0.06 0.06  1e3 0.0001]; 
ub=[0.6  4 0.8 0.8 0.3  0.3  0.1  0.1  1e7 10];
para=fmincon(@seir.loss,para,[],[],[],[],lb,ub); %%% 理论值比真实值最近
% para=fmincon(@seir.loss,para,[],[],[],[],[],[]);

Pi=para(3);
N0=para(9);
%%%%%%%%%%%%%%%%%%%% plot figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S0=N0-E0-I0-A0-R0;
ini=[S0 E0 I0 A0 R0];
[t,p]=ode45(@seir.dynamS,tspan,ini,[],para);



figure(1)
hold on;
scatter(t0,Ipresent(:,1),'.','r');
plot(t,p(:,3),'r');

%scatter(t0,Ipresent(:,2),'.','b');
%plot(t,p(:,5),'b');

%scatter(t0,Ipresent(:,3),'.','g');
%plot(t,p(:,4),'g');
 
xlabel('Time/day');
ylabel('Number of infections');
title('Fitting in Beijing')
%legend('Existing real infections','Fitting infections');
%figure(4)
hold on;
scatter(t,num(:,2),'.','b')
plot(t,p(:,4),'b');
 
%xlabel('Time/day');
%ylabel('Number of infections');
legend('Existing symptomatic infections','Fitting symptomatic infections','Existing asymptomatic infections','Fitting asymptomatic infections');



%%%%%%%%%%%%%% show the result %%%%%%%%%%%%%%%%%
c1=para(1)
c2=para(2)
Pi=para(3)
beta=para(4)
alpha1=para(5)
alpha2=para(6)
gama1=para(7)
gama2=para(8)
N=para(9)
decays=para(10)

%%%%%%%%%%%%%%%%%% show the beta(t) %%%%%%%%%%%%
figure(2)
hold on;
plot(t,c1*beta*exp(-decays*t)*Pi,'r');
plot(t,c2*beta*exp(-decays*t)*(1-Pi),'b');
xlabel('Time/day');
title('Beijing')
%ylabel('beta');
ylabel('Dynamic transmission parameters');%动态传播系数
legend('Dynamic transmission parameters of symptomatic infections','Dynamic transmission parameters of asymptomatic infections');
%figure(6)
%plot(t,c2*beta*exp(-decays*t)*(1-Pi));
%xlabel('时间/天');
%ylabel('beta');
%ylabel('无症状感染者传播系数');
%%%%%%%%%%%%%%%%%% show the E1 and E2 %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% show the E1 and E2 %%%%%%%%%%%%
figure(3)
hold on;
plot(t,c1*beta*exp(-decays*t)*Pi.*p(:,1).*p(:,2)/N,'r');
plot(t,c2*beta*exp(-decays*t)*(1-Pi).*p(:,1).*p(:,2)/N,'b');
title('Beijing')
xlabel('Time/day');
ylabel('Number of people infected / day');
legend('Number of people infected by symptomatic patients','Number of people infected by asymptomatic patients');
