function x0 = fitErrorFunc(X,Y,baseline)


nx = length(X);
Xrange = X(end)-X(1);
dX_mean = Xrange/nx;
ampl_guess = max(Y)-min(Y);
x0_guess = (X(1)+X(end))/2;
sigma_guess = Xrange/5;
if nargin == 3
    f_function = @(p,x) p(1)*0.5*(1+erf((x-p(2))/(p(3)*sqrt(2)))) + baseline;
else
    f_function = @(p,x) p(1)*0.5*(1+erf((x-p(2))/(p(3)*sqrt(2)))) + p(4);
end
baseline_guess = min(Y);
p_fit = lsqcurvefit(f_function,[ampl_guess,x0_guess,sigma_guess,baseline_guess],X,Y);
f_out = f_function(p_fit,X);
x0 = p_fit(2);
sigma_rise = p_fit(3);
is_x0_inRange = (x0 > X(1)) & (x0 < X(end));
if ~is_x0_inRange
    disp('Warning: fitted t0 is outside data range');
end
% xnew = X - x0;

% Plot result of t=0 fit
figure();
plot(X,Y,X,f_out);
legend('data','fit');
fit_results_str = {['t0 = ',num2str(round(x0)),' fs'],['\sigma = ',num2str(round(sigma_rise)), ' fs']};
dim = [0.7, 0.1, 0.2, 0.2];
annotation('textbox',dim,'String',fit_results_str,'FitBoxToText','on');
end