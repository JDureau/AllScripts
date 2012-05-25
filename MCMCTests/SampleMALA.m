function xstar = SampleMALA(x,Parameters)




Epsil = Parameters.Epsil;

epsilon = 10^(-10);
fx = log(max(eps,Parameters.f(x,Parameters)));
Grad = zeros(length(x),1);
for i = 1:Parameters.Dim
    xpdx = x;
    xpdx(i) = xpdx(i)+epsilon;
    fxpdx = log(max(eps,Parameters.f(xpdx,Parameters)));
    Grad(i,1) = (fxpdx-fx)/epsilon;
end

xstar = mvnrnd(x'+Epsil^2/2*Parameters.ScalingCov*Grad,squeeze(Epsil^2*Parameters.ScalingCov));

if isnan(xstar)
    disp('Pb Sample MALA')
end

% disp(mean(Grad))
% pause(0.01)


% plot(0,0,'+k')
% hold on
% plot(x(1),x(2),'ob')
% plot(xstar(1),xstar(2),'og')
% hold off
% xlim([-5 5])
% ylim([-1 1])
% pause()
