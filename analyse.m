function analyse(filename,dimensions);


filedump = csvread(filename);
configdata = filedump(1:19,6);
ndata = filedump(:,1:4);
gammadata = filedump(:,5);
timesteps = configdata(3);
N = configdata(2);
vol=N^dimensions;
t = (1:timesteps)';
tt = (0:timesteps)';


%--------------------------------------------------------------------------
% CURVE FITTING UNNECCESSARY WHEN LARGE NUMBERS OF SIM. RESULTS ARE
% AVERAGED TOGETHER TO PRODUCE NDATA & GAMMADATA
%
%POLYNOMIAL CURVE FITTING OF GAMMADATA
%
%figure(1);
%hold on;
%plot(gammadata);
%cubicfit=polyfit(t,gammadata,4);
%plot(t,polyval(cubicfit,t),'r')
%gammadata=polyval(cubicfit,t);
%
%NON LINEAR CURVE FITTING OF GAMMADATA
%
%figure(1);
%lb=-inf*ones(4,1);
%ub=inf*ones(4,1);
%params=[1,1000,0.004,1000];
%gaminterp = params(1)*t-params(2)*exp(-params(3)*t)+params(4);
%plot(gammadata);
%hold on;
%plot(t,gaminterp,'g');
options = optimset('LargeScale','off','TolX',1e-10,'TolFun',1e-10);
%[params,resnorm] = lsqcurvefit(@nonlinfun,[1,1000,0.004,1000],t,gammadata,lb,ub);
%gaminterp = params(1)*t-params(2)*exp(-params(3)*t)+params(4);
%plot(t,gaminterp,'r');
%
%NONLINEAR CURVE FITTING OF NDATA
%
%ninterp = zeros(timesteps,4);
%
%for i=2:3
%    if i == 3, params=[1,1000,0.004,1000]; end;
%    if i == 2, params=[-0.3,-0.7*configdata(8)*configdata(2)^2,0.003,0.3*configdata(8)*configdata(2)^2]; end;
%    ninterp(:,i) = params(1)*t-params(2)*exp(-params(3)*t)+params(4);
%    figure(2);
%    plot(ndata(:,i));
%    hold on;
%    plot(t,ninterp(:,i),'g');
%    [params,resnorm] = lsqcurvefit(@nonlinfun,params,t,ndata(:,i),lb,ub)
%    ninterp(:,i) = params(1)*t-params(2)*exp(-params(3)*t)+params(4);
%    plot(t,ninterp(:,i),'r');    
%    for j=2:timesteps-1
%        nderiv(j,i) = (ndata(j+1,i)-ndata(j-1,i))/2;
%    end
%    nderiv(1,i) = ndata(2,i) - ndata(1,i);
%    nderiv(timesteps,i) = ndata(timesteps,i) - ndata(timesteps-1,i);
%end;
%
%POLYNOMIAL CURVE FITTING OF NDATA
%
%for i=2:3
%    figure(2);
%    hold on;
%    plot(ndata(:,i));
%    cubicfit=polyfit(t,ndata(:,i),4);
%    plot(t,polyval(cubicfit,t),'r')
%    ndata(:,i)=polyval(cubicfit,t);
%end;
%
% ------------------------------------------------------------------------


% NUMERICALLY DIFFERENTIATE GAMMADATA & NDATA

gamderiv = zeros(timesteps,1);
for j=2:timesteps-1
    gamderiv(j) = (gammadata(j+1)-gammadata(j-1))/2;
end
gamderiv(1) = gammadata(2)-gammadata(1);
gamderiv(timesteps) = gammadata(timesteps) - gammadata(timesteps-1);

nderiv = zeros(timesteps,4);
for j=2:timesteps-1
        nderiv(j,:) = (ndata(j+1,:)-ndata(j-1,:))/2;
end
nderiv(1,:) = ndata(2,:) - ndata(1,:);
nderiv(timesteps,:) = ndata(timesteps,:) - ndata(timesteps-1,:);


% CALCULATE kX AT EACH TIMESTEP AND PLOT

k1 = vol * gamderiv ./ (ndata(:,1) .* ndata(:,2));
kNEG1 = (gamderiv + nderiv(:,2)) ./ ndata(:,4);
k2 = nderiv(:,3) ./ ndata(:,4);

figure(4);
subplot(2,1,1)
loglog(k1);
hold on;
subplot(2,1,2)
plot(kNEG1);
hold on;
plot(k2);


% FIT k(t)=a*t^(-h) MODEL TO k1 DATA USING LINEAR LEAST SQUARES
% PLOT AND WRITE h VALUE TO FILE.

figure(5);
hold on;
logk1=log(k1);
logt=log(t);
plot(logt(1:timesteps),logk1(1:timesteps));
C = ones(timesteps,2);
C(:,2) = logt(1:timesteps);
params = C\logk1(1:timesteps)
plot(logt(1:timesteps),params(2)*logt(1:timesteps)+params(1),'g');
configdata(12) = exp(params(1));
configdata(13) = -params(2);
h = configdata(13)
k0 = configdata(12)

% -------------------------------------------------------------------------
% FIT k(t)=a*(t-b)^(-h) MODEL TO k1 DATA USING NON-LINEAR LEAST SQUARES
% 
%
figure(3);
params=[3.75,2,0.5];
k1interp = params(1)*((t+params(2)).^(-params(3)));
%loglog(t(1:600),k1interp(1:600),'g');
%lb=[1,-100,-1];
%ub=[100,-0.1,1];
lb=[0,0,0];
ub=[Inf,Inf,1];
[params,resnorm] = lsqcurvefit(@k1fit,params,t,k1,lb,ub,options)
k1interp = params(1)*((tt+params(2)).^(-params(3)));
plot(tt,k1interp,'r');
hold on;
plot(t,k1,'b');
configdata(14:16)=real(params);
configdata(17)=configdata(14)/(configdata(15)^configdata(16))
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% FIT k(t)=a*t^(-h) MODEL TO k1 DATA USING NON-LINEAR LEAST SQUARES
% 
%
figure(3);
params=[3.75,0.5];
k1interp = params(1)*(t.^(-params(2)));
%loglog(t(1:600),k1interp(1:600),'g-');
%lb=[1,-100,-1];
%ub=[100,-0.1,1];
lb=-Inf;
ub=Inf;
[params,resnorm] = lsqcurvefit(@k1fitsimple,params,t(1:timesteps),k1(1:timesteps),lb,ub,options)
k1interp = params(1)*(t.^(-params(2)));
plot(t(1:timesteps),k1interp(1:timesteps),'g');
configdata(18:19)=real(params);
%
% -------------------------------------------------------------------------

% WRITE RESULTS TO FILE

filedump(1:19,6)=configdata;
filedump(:,7)=k1;
filedump(:,8)=kNEG1;
filedump(:,9)=k2;

csvwrite(filename,filedump);



