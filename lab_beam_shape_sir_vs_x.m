% построение зависимости SIR от расстояния для сценария: 1) сравнения
% управления формой луча с максимумом и максимумом и нулем ДН при RMSE=0 м;
% 2) сравнения управления формой луча с максимумом и максимумом и нулем ДН 
% при RMSE=5 м; 3) адаптивного управления формой/шириной ДН.
clear all; close all; clc; 
% выбор типа антенной решетки
% 1 - planar or uniform rectangural antenna array (URA), планарная АР
% 2 - uniform linear antenna array (ULA), линейная АР
% 3 - uniform circular antenna array (UCA); круговая АР: antPattCntrl=0,1,2
antType = 1;
Nel = 20;               % число АЭ в одном измерении
stdCoordsArr = [0 10];  % вектор RMSE (СКО) оценки координат UE, м
% вектор моделей управления формой и шириной луча при ДО:
%   0 - управление формой луча: максимумом ДН 
%   1 - управление формой луча: максимумом и нулем ДН 
%   2 - адаптивное управление формой/шириной ДНА (не зависит от RMSE)
%   3 - управление шириной луча 
antPattCntrlArr = [0, 1];

ZZa = [];
for s=1:length(stdCoordsArr)
    ZZ = [];
    stdCoords = stdCoordsArr(s);
    for aa=antPattCntrlArr
        antPattCntrl = aa;
        if aa ~= 2
            rng('default');
        else
            rng(s);
        end
        [X,Y,Z] = lab_beam_shape_fcn(antType,Nel,stdCoords,antPattCntrl);
        ZZ = [ZZ; Z];
    end
    ZZa = [ZZa; Z];
    figure; hold on;
    for ii=1:length(antPattCntrlArr)
        plot(X, ZZ(ii,:),'-', 'linewidth',2);
    end
    grid on; xlabel('x, m'); ylabel('Мгновенное SIR, дБ'); axis('tight');
    legend('управление максимумом ДН',...
           'управление максимумом и нулем ДН');
    title(sprintf("URA %i×%i, RMSE = %i м",Nel, Nel, stdCoords));
end
figure; hold on; axis('tight');
for ii=1:size(ZZa,1)
    plot(X, ZZa(ii,:),'-','linewidth',2);
end
grid on; xlabel('x, м'); ylabel('SIR, дБ');
title(sprintf("управление максимумом и нулем ДН, URA %i×%i",Nel, Nel));
legend(sprintf("RMSE = %i м",stdCoordsArr(1)),...
       sprintf("RMSE = %i м",stdCoordsArr(2)));
axis('tight');