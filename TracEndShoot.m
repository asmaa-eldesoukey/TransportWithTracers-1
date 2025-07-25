function Res = EndShoot(PI,PhiF,T,St,dSt,n)
y0 = [1;0;0;1;PI];
[~,y] = ode45(@(t,y) TracODE(t,y,St,dSt,n),T,y0);
PhiFs = reshape(y(end,1:n*n),[n,n]);

Res = reshape(PhiFs-PhiF,[n^2,1]);
end

U = [cos(a) sin(a);
    -sin(a) cos(a)]
PhiF = sqrtm(SigmaF)*U;


cond = 1;
N = 1;

options = optimoptions("fsolve","MaxFunctionEvaluations",1e4,...
    "Algorithm",'levenberg-marquardt');

while 1
if (N<10000 && cond>=1E-3)
   P0 = randn(n*n,1);
   F = @(P) EndShoot(P,PhiF,T,Sigma,dotSigma,n);
   [P,fval] = fsolve(F,P0,options);
   cond = norm(fval);
   N = N+1;
else
     display(P);
     display(fval);
     display(N)
    break
    
end
end

y0 = [1;0;0;1;P];
sole = ode45(@(t,y) TracODE(t,y,Sigma,dotSigma,n),T,y0);

figure,plot(sole.x,sole.y(1,:),sole.x,sole.y(2,:),...
    sole.x,sole.y(3,:),sole.x,sole.y(4,:),LineWidth=2)
ax = gca;
ax.DataAspectRatio=[0.1,1,1];
ax.FontSize = 15;
legend('$\Phi_t^\star(1,1)$','$\Phi_t^\star(2,1)$','$\Phi_t^\star(1,2)$','$\Phi_t^\star(2,2)$',...
                    'FontSize',14,'interpreter','latex','Location','northeastoutside');
xlabel('Time \ t', 'FontSize',19 ,'interpreter','latex');

Phit_mat = zeros(n,n,size(sole.y,2));
sigma_out = zeros(n,n,size(sole.y,2));
pos_out = zeros(2,length(sole.x));
posi = PosI;

for k = 1: size(Phit_mat,3)
    Phit_mat(:,:,k) = [sole.y(1,k),sole.y(3,k);
                       sole.y(2,k),sole.y(4,k)];
    sigma_out(:,:,k) =  Phit_mat(:,:,k)*Phit_mat(:,:,k)';
     pos_out(:,k) = Phit_mat(:,:,k)*posi;
end


angle_out = acosd(dot(pos_out(:,1),pos_out(:,end))./(norm(pos_out(:,1))*norm(pos_out(:,end))));
err_dual = SigmaF-Phit_mat(:,:,end)*Phit_mat(:,:,end)';


resol = 100;
xo = zeros(2,resol);
yo = zeros(2,resol);
curv = linspace(0,2.*pi,resol);
v = linspace(0,1);
Theta = linspace(0,2*pi,1e2);

sigma__out_(:,:,1) = sigma_out(:,:,1);
sigma__out_(:,:,2) = sigma_out(:,:,end);



L1 =  reshape(0.5.*(sigma__out_(1,1,:)+sigma__out_(2,2,:))+...
           sqrt((0.5.*sigma__out_(1,1,:)-0.5.*sigma__out_(2,2,:)).^2+...
           sigma__out_(1,2,:).^2),1,[]);
L2 = reshape(0.5.*(sigma__out_(1,1,:)+sigma__out_(2,2,:))-...
           sqrt((0.5.*sigma__out_(1,1,:)-0.5.*sigma__out_(2,2,:)).^2+...
           sigma__out_(1,2,:).^2),1,[]);
theta = zeros(1,Nt);


for tt = 1:2
    if (sigma__out_(1,2,tt)==0 && sigma__out_(1,1,tt)>=sigma__out_(2,2,tt))
            theta(tt) = 0;
    elseif (sigma__out_(1,2,tt)==0 && sigma__out_(1,1,tt)<sigma_out(2,2,tt))
         theta(tt) = pi/2;
    else
        theta(tt) = atan2(L1(tt)-sigma__out_(1,1,tt),sigma__out_(1,2,tt));
    end
end




for k = 1:resol
    for i = 1:2
    xo(i,k) = sqrt(L1(i)).*cos(theta(i)).*cos(curv(k))-...
        sqrt(L2(i)).*sin(theta(i)).*sin(curv(k));
    yo(i,k) = sqrt(L1(i)).*sin(theta(i)).*cos(curv(k))+...
        sqrt(L2(i)).*cos(theta(i)).*sin(curv(k));
    end
end



Xo = zeros(length(sole.x),length(Theta));
Yo = zeros(length(sole.x),length(Theta));
Zo = zeros(length(sole.x),length(Theta));

for i =1:length(T)
    for j = 1:length(Theta)
        s = Theta(j);
        Xo(i,j) = T(i);
        Yo(i,j) = [1,0]*(sqrtm(Sigma(T(i))))*[cos(s);sin(s)];
        Zo(i,j) = [0,1]*(sqrtm(Sigma(T(i))))*[cos(s);sin(s)];
    end
end


figure,surf(Xo,Yo,Zo,[0 1 0 2*pi],'LineWidth',0.1,'FaceColor','k',...
    'FaceAlpha',0.1,'MeshStyle','none'), hold on;
plot3(0.*ones(resol,1),xo(1,:),yo(1,:),'k',LineWidth=2); hold on;
fill3(1.*ones(resol,1),xo(end,:),yo(end,:),'k',FaceAlpha=0.125,EdgeColor='none',LineWidth=2);hold on
% plot3(1.*ones(50,1),xo(end,1:50),yo(end,1:50),'k',LineWidth=2); hold on;
% plot3(1.*ones(50,1),xo(end,51:100),yo(end,51:100),'k--',LineWidth=2);
plot3(1.*ones(100,1),xo(end,:),yo(end,:),'k',LineWidth=2)
hold on
   plot3(sole.x,pos_out(1,:),pos_out(2,:),'Color',[0.4660 0.6740 0.1880],LineWidth=3);
   hold on;
plot3(0.*ones(resol,1),linspace(0,pos_out(1,1),resol),...
    linspace(0,pos_out(2,1),resol),':','Color', [0.41 0.41 0.41],LineWidth=2); hold on;
plot3(ones(resol,1,1),linspace(0,pos_out(1,end),resol),....
    linspace(0,pos_out(2,end),resol),':','Color', [0.41 0.41 0.41],LineWidth=2); hold on;
plot3(linspace(0,1,resol),zeros(resol,1),zeros(resol,1),...
    ':','Color', [0.41 0.41 0.41],LineWidth=2);hold on;
ax = gca;
ax.DataAspectRatio=[0.125,1,1];
ax.FontSize = 15;
xlabel('Time \ t', 'FontSize',19 ,'interpreter','latex');
% view(-90,0)

