function params=fitSymTilt111(stgb111,params)
ergb=params(1);
ang=pi/3*180/pi;
C=stgb111(abs(stgb111(:,1)-ang)==min(abs(stgb111(:,1)-ang)),2);
a=[params(40) params(39)];
a(1)=a(1)*ergb;

da=0.1;
theta=stgb111(:,1)';
E=stgb111(:,2)';

[theta, i]=sort(theta);
E=E(i);

A=stepParameters(C,a,da,theta,E);
A=[A C];

dx=(max(theta)-min(theta))/200;
x=min(theta):dx:max(theta);
y=equation(x,A);
h=figure(1);
plot(theta,E,'.',x,y,'-')

A(1)=A(1)/ergb;
A(3)=A(3)/ergb;

title(strcat('UO2 111 Symmetric Tilt    eRGB: ',num2str(ergb),' Jm^{-2}'))
text(20,1.2,strcat('Max Energy: ', num2str(A(1))))
text(20,1.1,strcat('Max Angle: ', num2str(A(2))))
text(20,1.0,strcat('Min Energy: ', num2str(A(3))))
xlabel('Angle')
ylabel('Energy')
saveas(h,'./Results/111SymTilt.jpg')

params(39)=A(2);
params(40)=A(1);
params(41)=A(3);

end

function a=stepParameters(C,a,da,x,y)
tol=0.000000001;
it=0;
while da>tol
    it=it+1;
    g=equation(x,[a C]);
    now=leastSquare(g,y);
    for i=1:length(a)
        a(i)=a(i)+da;
        g=equation(x,[a C]);
        X2(i,1)=leastSquare(g,y);
        a(i)=a(i)-da;
    end
    for i=1:length(a)
        a(i)=a(i)-da;
        g=equation(x,[a C]);
        X2(i,2)=leastSquare(g,y);
        a(i)=a(i)+da;
    end
    
    [val ind]=min(X2(:));
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
    
    if mod(it,5)==0
        dx=(max(x)-min(x))/200;
        px=min(x):dx:max(x);
        py=equation(px,[a C]);
        figure(1)
        plot(x,y,'.',px,py,'-')
        title('111 Symmetric Tilt')
        xlabel('Misorientation Angle')
        ylabel('Energy')
        drawnow
    end
end


end

function E=equation(x,a)
x=x*pi/180;

x(x > pi/3) = 2*pi/3 - x(x>pi/3);

en1=0;
theta1=0;
theta2=pi/3;

ksim=a(2);
enmax=a(1);
enmin=a(3);

a1=0.5;
a2=0.5;

cx1=x>=theta1;
cx2=x<ksim;
cx=cx1&cx2;
dtheta=ksim-theta1;
theta=x(cx);
thetan=(theta-theta1)./dtheta*pi/2;
E1=en1+(enmax-en1).*(sin(thetan)-a1.*(sin(thetan).*log(sin(thetan))));
E1(thetan==0)=en1;

cx1=x>=ksim;
cx2=x<=theta2;
cx=cx1&cx2;
dtheta=ksim-theta2;
theta=x(cx);
thetan=(theta-theta2)./dtheta*pi/2;
E2=enmin+(enmax-enmin).*(sin(thetan)-a2.*(sin(thetan).*log(sin(thetan))));
E2(thetan==0)=enmin;

E=[E1 E2];
end

function X2=leastSquare(g,y)
X2=sum((y-g).^2);
end