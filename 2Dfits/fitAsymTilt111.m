function params=fitAsymTilt111(atgb111,params)
    ksi=atgb111(:,1)';
    eta=atgb111(:,2)';
    E=atgb111(:,3)';
    
    [ksi, i]=sort(ksi);
    eta=eta(i);
    E=E(i);
    
    da=0.1;
    a=[params(42) params(43)];
    
    A=stepParameters(a,da,ksi,eta,E,params);
    
    params(42)=A(1);
    params(43)=A(2);
end

function a=stepParameters(a,da,ksi,eta,y,params)
tol=0.0000001;
it=0;
while da>tol
    it=it+1;
    g=equation(ksi,eta,params,a);
    now=leastSquare(g,y);
    for i=1:length(a)
        a(i)=a(i)+da;
        g=equation(ksi,eta,params,a);
        X2(i,1)=leastSquare(g,y);
        a(i)=a(i)-da;
    end
    for i=1:length(a)
        a(i)=a(i)-da;
        g=equation(ksi,eta,params,a);
        X2(i,2)=leastSquare(g,y);
        a(i)=a(i)+da;
    end
    
    [val ind]=min(X2(:));
    
    if mod(it,10)==0
        k=0:5:60;
        e=0:5:60;
        [K,E]=meshgrid(k,e);
        for i=1:length(k)
            for j=1:length(e)
                p(i,j)=equation(k(i),e(j),params,a);
            end
        end
        figure(1);
        surf(E,K,p)
        title('111 Asymmetric Tilt')
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
k=0:1:120;
e=0:1:120;
[K,E]=meshgrid(k,e);
for i=1:length(k)
    for j=1:length(e)
        p(i,j)=equation(k(i),e(j),params,a);
    end
end
h=figure(1);
surf(E,K,p)
title(strcat('111 Asymmetric Tilt   Symmetry Point:',num2str(a(1)),'   Eta Scale:',num2str(a(2))))
xlabel('Ksi')
ylabel('Eta')
zlabel('Energy')
hold on;
scatter3(ksi,eta,y,40,'m','filled');
drawnow;
hold off;
pause;
saveas(h,'./Results/111AsymTilt.jpg')
end

function en=equation(ksi,eta,params,a)
    % en = atgbs111(eta,ksi,pars)
%
% This function is a fit to the energies of all 111-tilt boundaries

% There's an additional symmetry in 111 atgbs that doesn't exist in 100 or
% 110 atgbs.  This is because a rotation about [111] equal to half the period
% (i.e. 60 degrees) is equivalent to a mirror reflection in the (111)
% plane.  Both are Sigma3 operations.  The same is not true of the
% 45-degree [100] or the 90-degree [110] rotation.
% The following two lines account for this extra symmetry.
    ksi=ksi*pi/180;
    eta=eta*pi/180;
    
    ksi(ksi > pi/3) = 2*pi/3 - ksi(ksi>pi/3);
    eta(eta > pi/3) = 2*pi/3 - eta(eta>pi/3);

    % Below the following value of ksi, we ignore the eta dependence.  This is
    % because there's little evidence that it actually varies.  Above this
    % value, we interpolate on an rsw function that follows the Sigma3 line,
    % which is also a line of symmetry for the function.
    ksim  = params(39); % 111 atgb ksi break

    enmax = params(40)*params(1); % Energy at the peak (ksi == ksim)
    enmin = params(41)*params(1); % energy at the minimum (sigma3, eta == 0)
    encnt = a(1)*params(1); % energy at the symmetry point (sigma3, eta == pi/3)

    a1    = 0.5;
    a2    = 0.5;

    etascale = a(2); % eta scaling parameter for 111 atgb rsw function on Sigma3 line
        % This rsw function is unusual in that the change in shape of the
        % function is much better captured by changing the angular scale rather
        % than changing the dimensionless shape factor.

    en = zeros(size(ksi));

    select = (ksi <= ksim);
    en(select)  = enmax*rsw(ksi(select),0,ksim,a1) ;

    % chi is the shape of the function along the sigma3 line.
    chi = enmin + (encnt-enmin)*rsw(eta(~select),0,pi/(2*etascale),0.5);
    en(~select)  = chi   + (enmax - chi).*rsw(ksi(~select),pi/3,ksim,a2) ;
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