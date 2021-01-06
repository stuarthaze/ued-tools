function checkFittingDiagnostics(FIT,options)

if options.saveFigs
    saveDirectoryName = options.saveDirectoryName;
    if ~isdir(['.\',saveDirectoryName])
        mkdir(saveDirectoryName);
    end
end
figure();
plot(FIT.exitFlags); title('Optimizer exitflags - Converged if >= 1');

Residuals2rd_initial = FIT.initialResiduals.^2;
Residuals2rd_refined = FIT.Residuals.^2;
figure();
subplot(1,2,1); imagesc(Residuals2rd_initial); colorbar(); title('Residuals^2 Initial');
subplot(1,2,2); imagesc(Residuals2rd_refined); colorbar(); title('Residuals^2 Refined');
if options.saveFigs
    savefig(['.\',saveDirectoryName,'\Residuals_2d']);
end

% Residuals
RMS_Residual_vs_T_initial = sqrt(mean(FIT.initialResiduals.^2,2));
RMS_Residual_vs_T_refined = sqrt(mean(FIT.Residuals.^2,2));
figure(); plot([RMS_Residual_vs_T_initial,RMS_Residual_vs_T_refined]); title('RMS residuals');
legend('Initial','Refined');
if options.saveFigs
    savefig(['.\',saveDirectoryName,'\Residuals_RMS_1d']);
end

% Restraint function values
figure(); plot(sum(FIT.RestraintVals,2)); title('Restraint function sum');
if options.saveFigs
    savefig(['.\',saveDirectoryName,'\Restaint_sum']);
end
%
% Uiso
if ~isempty(FIT.Ufree_T)
    figure(); plot(FIT.Ufree_T'); title('Uiso fitted'); 
    if options.saveFigs
        savefig(['.\',saveDirectoryName,'\Uiso']);
    end
end
% DWF
figure();
plot(FIT.DW_es);
% hold on;
% plot(FIT.DW_gsHot);
title('Debye-Waller exponent <d(u^2)>');
% legend('ES','hot GS');
if options.saveFigs
    savefig(['.\',saveDirectoryName,'\Debye-Waller Factor']);
end

% Position changes
meanPositionChanges = mean(FIT.PositionChanges,1);
maxPositionChanges = max(FIT.PositionChanges,[],1);

figure(); plot([meanPositionChanges', maxPositionChanges']); 
title('dX'); legend('mean', 'max');
if options.saveFigs
    savefig(['.\',saveDirectoryName,'\Position_changes_mean_max']);
end
end
