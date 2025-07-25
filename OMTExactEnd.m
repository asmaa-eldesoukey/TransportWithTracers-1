clc;
clear;
close all;

Nt = 400;
T = linspace(0,1,Nt);
dt = T(2)-T(1);
n = 2;
g = 0;
PosI = [-cos(g);sin(g)];
a = 67.5*pi/180;
Nmax = 10000;


SigmaF = [2,sqrt(2);sqrt(2),2];
PhioF = sqrtm(SigmaF);


Phio = @(t) (1-t).*eye(2) + t.*PhioF;
Sigma = @(t) Phio(t)*Phio(t)';
A = @(t)(PhioF-eye(2))*pinv(eye(2)+t.*(PhioF-eye(2))); 
dotSigma = @(t) A(t)*Sigma(t)+Sigma(t)*A(t)';



PhioEval = zeros(n,n,Nt);


for i = 1:Nt
    PhioEval(:,:,i) =Phio(T(i));
end


figure,plot(T,squeeze(PhioEval(1,1,:)),T,squeeze(PhioEval(1,2,:)),....
    T,squeeze(PhioEval(2,1,:)),T,squeeze(PhioEval(2,2,:)),LineWidth=2)
legend('$\Phi_{11}^{\rm o}$','$\Phi_{12}^{\rm o}$',...
    '$\Phi_{21}^{\rm o}$','$\Phi_{22}^{\rm o}$','FontSize',19,'interpreter','latex');

