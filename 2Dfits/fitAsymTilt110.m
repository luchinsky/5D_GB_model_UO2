function params=fitAsymTilt110(atgb110,params)
    ksi=atgb110(:,1)';
    eta=atgb110(:,2)';
    E=atgb110(:,3)';
    
    [ksi, i]=sort(ksi);
    eta=eta(i);
    E=E(i);
    
    eta0=symTilt110(ksi,params);
    eta180=symTilt110(180-ksi,params);
    
    da=0.1;
    a=params(26);
    
    A=stepParameters(a,da,eta0,eta180,ksi,eta,E,params);
    
    params(26)=A;
end

function a=stepParameters(a,da,eta0,eta180,ksi,eta,y,params)
tol=0.0000001;
it=0;
while da>tol
    it=it+1;
    g=equation(eta,eta0,eta180,a);
    now=leastSquare(g,y);
    for i=1:length(a)
        a(i)=a(i)+da;
%         if a(i)>1
%             a(i)=1;         %artificially set maximum to 1
%             g=equation(eta,eta0,eta180,a);
%             X2(i,1)=leastSquare(g,y);
%         else
            g=equation(eta,eta0,eta180,a);
            X2(i,1)=leastSquare(g,y);
            a(i)=a(i)-da;
%         end
    end
    for i=1:length(a)
        a(i)=a(i)-da;
        g=equation(eta,eta0,eta180,a);
        X2(i,2)=leastSquare(g,y);
        a(i)=a(i)+da;
    end
    
    [val ind]=min(X2(:));
    
    if mod(it,10)==0
        k=0:5:180;
        e=0:5:180;
        [K,E]=meshgrid(k,e);
        peta0=symTilt110(k,params);
        peta180=symTilt110(180-k,params);
        for i=1:length(k)
            for j=1:length(e)
                p(i,j)=equation(e(j),peta0(i),peta180(i),a);
            end
        end
        figure(1);
        surf(E,K,p)
        title('110 Asymmetric Tilt')
        xlabel('Ksi')
        ylabel('Eta')
        zlabel('Energy')
        drawnow;
    end
    
    col=1;
    if ind>length(a)
        ind=ind-length(a);
        col=2;
    end
    if now<=min(X2(:))
        da=da/2;
    elseif col==1
        a(ind)=a(ind)+da;
    elseif col==2
        a(ind)=a(ind)-da;
    end
    
end
k=0:2:180;
e=0:2:180;
[K,E]=meshgrid(k,e);
peta0=symTilt110(k,params);
peta180=symTilt110(180-k,params);
for i=1:length(k)
    for j=1:length(e)
        p(i,j)=equation(e(j),peta0(i),peta180(i),a);
    end
end
h=figure(1);
surf(E,K,p)
title(strcat('110 Asymmetric Tilt   Shaping Factor:',num2str(a(1))))
xlabel('Ksi')
ylabel('Eta')
zlabel('Energy')
hold on;
scatter3(ksi,eta,y,40,'m','filled');
drawnow;
hold off;
saveas(h,'./Results/110AsymTilt.jpg')
pause
end

function E=symTilt110(x,params)
x=x*pi/180;

en1=0;
theta1=0;
theta3=acos(1/3);
theta5=acos(-7/11);
theta7=pi;
en7=en1;


en2=params(27)*params(1);
en3=params(28)*params(1);
en4=params(29)*params(1);
en5=params(30)*params(1);
en6=params(31)*params(1);
theta2=params(32);
theta4=params(33);
theta6=params(34);

a1=0.5;
a2=0.5;
a3=0.5;
a4=0.5;
a5=0.5;
a6=0.5;

thetan=zeros(1,length(x));
E1=zeros(1,length(x));
E2=zeros(1,length(x));
E3=zeros(1,length(x));
E4=zeros(1,length(x));
E5=zeros(1,length(x));
E6=zeros(1,length(x));

cx1=x>=theta1;
cx2=x<theta2;
cx=cx1&cx2;
dtheta=theta2-theta1;
thetan(cx)=(x(cx)-theta1)./dtheta*pi/2;
E1(cx)=en1+(en2-en1).*(sin(thetan(cx))-a1.*(sin(thetan(cx)).*log(sin(thetan(cx)))));
E1(isnan(E1))=en1;

cx1=x>=theta2;
cx2=x<theta3;
cx=cx1&cx2;
dtheta=theta2-theta3;
thetan(cx)=(x(cx)-theta3)./dtheta*pi/2;
E2(cx)=en3+(en2-en3).*(sin(thetan(cx))-a2.*(sin(thetan(cx)).*log(sin(thetan(cx)))));
E2(isnan(E2))=en3;

cx1=x>=theta3;
cx2=x<theta4;
cx=cx1&cx2;
dtheta=theta4-theta3;
thetan(cx)=(x(cx)-theta3)./dtheta*pi/2;
E3(cx)=en3+(en4-en3).*(sin(thetan(cx))-a3.*(sin(thetan(cx)).*log(sin(thetan(cx)))));
E3(isnan(E3))=en3;

cx1=x>=theta4;
cx2=x<theta5;
cx=cx1&cx2;
dtheta=theta4-theta5;
thetan(cx)=(x(cx)-theta5)./dtheta*pi/2;
E4(cx)=en5+(en4-en5).*(sin(thetan(cx))-a4.*(sin(thetan(cx)).*log(sin(thetan(cx)))));
E4(isnan(E4))=en5;

cx1=x>=theta5;
cx2=x<theta6;
cx=cx1&cx2;
dtheta=theta6-theta5;
thetan(cx)=(x(cx)-theta5)./dtheta*pi/2;
E5(cx)=en5+(en6-en5).*(sin(thetan(cx))-a5.*(sin(thetan(cx)).*log(sin(thetan(cx)))));
E5(isnan(E5))=en5;

cx1=x>=theta6;
cx2=x<theta7;
cx=cx1&cx2;
dtheta=theta6-theta7;
thetan(cx)=(x(cx)-theta7)./dtheta*pi/2;
E6(cx)=en7+(en6-en7).*(sin(thetan(cx))-a6.*(sin(thetan(cx)).*log(sin(thetan(cx)))));
E6(isnan(E6))=en7;

E=E1+E2+E3+E4+E5+E6;
%E=[E1 E2 E3 E4 E5 E6];
end

function en=equation(eta,en1,en2,a)
    eta=eta*pi/180;
    
    en = zeros(size(eta));

    % Power-law interpolation did not work well in this case.  Did an rsw
    % function instead.
    select = en1>=en2;

    en(select) = en2(select) + (en1(select)-en2(select)).*rsw(eta(select),pi,0,a);
    en(~select) = en1(~select) + (en2(~select)-en1(~select)).*rsw(eta(~select),0,pi,a);

end

function X2=leastSquare(g,y)
X2=sum((y-g).^2);
end

function en = rsw(theta,theta1,theta2,a)
% en = rsw(theta,theta1,theta2,a)
%
% This function computes the value of Read-Shockley-Wolf function at theta.
% The RSW function is normalized to be 1.0 at theta2 and 0.0 at theta1.
% 
% theta             angle at which to compute the function
% theta1            the starting angle of the interval
% theta2            the end angle of the interval
% a                 parameter defining the shape of the RSW function
%

    dtheta = theta2 - theta1  ;     % Interval of angles where defined
    theta = (theta-theta1)./dtheta*pi/2 ;    % Normalized angle
    % The rest is the RSW function evaluation
    sins = sin(theta) ;
    xlogx = zeros(size(sins));

    % Cut off at small sins to avoid 0*infinity problem.  The proper limit is 0.
    select = sins >= 0.000001;
    xlogx(select) = sins(select).*log(sins(select));

    en = sins - a*xlogx ;
end