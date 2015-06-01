function params=fitTwist100(twist100,params)
ergb=params(1);
a=[params(10) params(11)];
a(1)=a(1)*ergb;

da=0.1;
theta=twist100(:,1)';
E=twist100(:,2)';
[theta, i]=sort(theta);
E=E(i);

A=stepParameters(a,da,theta,E);

dx=(max(theta)-min(theta))/200;
x=min(theta):dx:max(theta);
y=equation(x,A);
h=figure(1);
plot(theta,E,'.',x,y,'-')

A(1)=A(1)/ergb;

title(strcat('UO2 100 Twist   eRGB: ',num2str(ergb),' Jm^{-2}'))
text(20,0.6,strcat('Max Energy: ', num2str(A(1))))
text(20,0.5,strcat('Shape Factor: ', num2str(A(2))))
xlabel('Angle')
ylabel('Energy')
saveas(h,'./Results/100Twist.jpg')

params(10)=A(1);
params(11)=A(2);

end

function a=stepParameters(a,da,x,y)
tol=0.000000001;
it=0;
while da>tol
    it=it+1;
    g=equation(x,a);
    now=leastSquare(g,y);
    for i=1:length(a)
        a(i)=a(i)+da;
        g=equation(x,a);
        X2(i,1)=leastSquare(g,y);
        a(i)=a(i)-da;
    end
    for i=1:length(a)
        a(i)=a(i)-da;
        g=equation(x,a);
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
        py=equation(px,a);
        figure(1)
        plot(x,y,'.',px,py,'-')
        title('100 Twist')
        xlabel('Misorientation Angle')
        ylabel('Energy')
        drawnow
    end
end


end

function E=equation(x,a)
x=x*pi/180;

en1=0;
theta1=0;
theta2=pi/4;

en2=a(1);
a=a(2);

a1=a*en2;

dtheta=theta2-theta1;
theta=x;
thetan=(theta-theta1)./dtheta*pi/2;
E1=en2.*sin(thetan)-a1.*(sin(thetan).*log(sin(thetan)));
E1(thetan==0)=en1;

E=E1;
end

function X2=leastSquare(g,y)
X2=sum((y-g).^2);
end