function params=fitWeight(rgb,params)
    da=0.1;
    a=[params(5) params(6) params(7)];
    
    A=stepParameters(a,da,rgb,params);
    A
    params(5)=A(1);
    params(6)=A(2);
    params(7)=A(3);
end

function a=stepParameters(a,da,rgb,params)
y=rgb(:,19);
tol=0.0001;
it=0;
maxits=30;
while da>tol
    it=it+1;
    g=equation(rgb,a,params);
    now=leastSquare(g,y);
    for i=1:length(a)
        a(i)=a(i)+da;
        g=equation(rgb,a,params);
        X2(i,1)=leastSquare(g,y);
        a(i)=a(i)-da;
    end
    for i=1:length(a)
        a(i)=a(i)-da;
        if a(i)<=0
            X2(i,2)=now*10;
        else
            g=equation(rgb,a,params);
            X2(i,2)=leastSquare(g,y);
        end
        a(i)=a(i)+da;
    end
    
    [val ind]=min(X2(:));
    
    pX2(it)=val;
    
    if mod(it,2)==0
        figure(1);
        plot(pX2);
        title(strcat('Weight Parameters   da:',num2str(da)))
        xlabel('Iterations')
        ylabel('Error')
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
    if it>maxits
        break
    end
end


end

function E=equation(rgb,a,params)
    [m,n]=size(rgb);
    params(5)=a(1);
    params(6)=a(2);
    params(7)=a(3);
    
    for i=1:m
        [P,Q]=convert2PQ(rgb(i,1:18));
        E(i)=GB5DOF(P,Q,params);
    end
    E=E';
end

function X2=leastSquare(g,y)
X2=sum((y-g).^2);
end