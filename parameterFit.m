function main()
    clear; clc;
    addpath(genpath(pwd));
    
    runPrediction=1;
    axis=[2 1 0];
    ksimax=180;
    phi=0;
    
    h=waitbar(0,'Reading in Files');
    getFitData()
    getTestData()
    getParameters()
    
    waitbar(0.1,h,'Fitting Twist Parameters')
    fitTwist()
    
    waitbar(0.2,h,'Fitting Symmetric Tilt Parameters')
    fitSymmetricTilt()
    
    waitbar(0.3,h,'Fitting Asymmetric Tilt Parameters')
    fitAsymmetricTilt();
    delete(h);
    
    fitRemaining()
    
    h=waitbar(0.9,'Testing Fit');
    testFit()
    
    if runPrediction==1
        h=waitbar(0.92,'Calculating Prediction');
        predict(axis,ksimax,phi)
    end
    
    waitbar(0.95,h,'Writing Parameters')
    writeParameters()
    
    waitbar(1,h,'Complete!')
    pause(0.5);
    delete(h);
end

function getFitData()
    global twist100;
    global twist110;
    global twist111;
    global stgb100;
    global stgb110;
    global stgb111;
    global atgb100;
    global atgb110;
    global atgb111;
    global rgb;
    
    fn100twist='./InputData/100Twist.xlsx';
    fn110twist='./InputData/110Twist.xlsx';
    fn111twist='./InputData/111Twist.xlsx';
    fn100stgb='./InputData/100SymmTilt.xlsx';
    fn110stgb='./InputData/110SymmTilt.xlsx';
    fn111stgb='./InputData/111SymmTilt.xlsx';
    fnatgb='./InputData/AsymmTilt.xlsx';
    fnrandom='./InputData/RandomGB.xlsx';
    
    twist100=xlsread(fn100twist);
    twist110=xlsread(fn110twist);
    twist111=xlsread(fn111twist);
    stgb100=xlsread(fn100stgb);
    stgb110=xlsread(fn110stgb);
    stgb111=xlsread(fn111stgb);
    asym=xlsread(fnatgb);
    rgb=xlsread(fnrandom);
    
    atgb=findKsiEta(asym);
    atgb100=[atgb(atgb(:,1)==100,2) atgb(atgb(:,1)==100,3) atgb(atgb(:,1)==100,5)];
    atgb110=[atgb(atgb(:,1)==110,2) atgb(atgb(:,1)==110,3) atgb(atgb(:,1)==110,5)];
    atgb111=[atgb(atgb(:,1)==111,2) atgb(atgb(:,1)==111,3) atgb(atgb(:,1)==111,5)];
    
    atgb100
    atgb110
    atgb111
end

function getTestData()
    global testData;
    
    fntest='./InputData/TestData.xlsx';
    
    testData=xlsread(fntest);
    
end

function getParameters()
    global params;
    
    fnParams='./Results/43Parameters.xlsx';
    
    params=xlsread(fnParams);
end

function writeParameters()
    global params;
    
    fnParams='./Results/43Parameters.xlsx';
    
    xlswrite(fnParams,params);
end

function testFit()
    global testData;
    global params;
    
    [m,n]=size(testData);
    X2=0;
    en=zeros(m,1);
    for i=1:m
        [P,Q]=convert2PQ(testData(i,:));
        en(i)=GB5DOF(P,Q,params);
        err(i)=(en(i)-testData(i,19))*100/testData(i,19);
    end
    
    H=figure(1);
    hist(err)
    title('Test Fit')
    xlabel('Percent Deviation from MD Values')
    ylabel('Number of Points')
    saveas(H,'./Results/Validation.jpg')
end

function fitTwist()
    global twist100;
    global twist110;
    global twist111;
    global params;
    
    params=fitTwist100(twist100,params);
    params=fitTwist110(twist110,params);
    params=fitTwist111(twist111,params);
end

function fitSymmetricTilt()
    global stgb100;
    global stgb110;
    global stgb111;
    global params;
    
    params=fitSymTilt100(stgb100,params);
    params=fitSymTilt110(stgb110,params);
    params=fitSymTilt111(stgb111,params);
end

function fitAsymmetricTilt()
    global atgb100;
    global atgb110;
    global atgb111;
    global params;
    params(12)
    params=fitAsymTilt100(atgb100,params);
    params(12)
    params(26)
    params=fitAsymTilt110(atgb110,params);
    params(26)
    params(42)
    params(43)
    params=fitAsymTilt111(atgb111,params);
    params(42)
    params(43)
end

function fitRemaining()
    global rgb;
    global params;
    
    h=waitbar(0.4,'Fitting Remaining Parameters');
    maxits=3;
    for i=1:maxits
        waitbar(0.3+0.15*i,h,'Fitting Mixing Parameters')
        params=fitMix(rgb,params);
        
        waitbar(0.4+0.15*i,h,'Fitting Weighing Parameters')
        params=fitWeight(rgb,params);
    end
    delete(h);
end

function predict(axis,ksimax,phi)
    global params;
    ang=[8.525386 16.95743 25.20876 31.23226 40.87863 48.18968 58.41186 73.39845 96.37937 131.8103 154.7912 167.2413 180];
    E=[1.588251 1.901526 2.014194 2.018589 2.016792 1.708756 2.051193 2.024049 2.102244 1.701356 1.993355 1.789079 1.657378];
    
    dk=ksimax/50;
    ksi=0:dk:ksimax;
    eta=0;
    
    for i=1:length(ksi)
        [P Q]=ksieta2PQ(axis, ksi(i), eta, phi);
        en(i)=GB5DOF(P,Q,params);
    end
    
    figure(1)
    plot(ksi,en,'-b',ang,E,'.b')
    title(strcat('Prediction for [',num2str(axis(1)),' ',num2str(axis(2)),' ',num2str(axis(3)),']   Phi:',num2str(phi)))
    xlabel('Angle')
    ylabel('Energy')
    drawnow;
end
