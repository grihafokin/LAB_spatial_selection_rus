clear all; close all; clc;
stdCoordsArr = [1, 3, 5, 10];  % СКО оценки координат по [x, y, z], м
NellArr = 4:4:20;              % число АЭ в одном измерении
% выбор типа антенной решетки:
% 1 - planar or uniform rectangural antenna array (URA), планарная АР
% 2 - uniform linear antenna array (ULA), линейная АР
% 3 - uniform circular antenna array (UCA); круговая АР: не поддерживается
antType = 1;

% вектор моделей управления формой и шириной луча при ДО:
%   0 - управление формой луча: максимумом ДН 
%   1 - управление формой луча: максимумом и нулем ДН 
%   2 - адаптивное управление формой/шириной ДНА (не зависит от RMSE)
%   3 - управление шириной луча 
%   адаптивное управление ДН не зависит от СКО и не рассматривается
antPattCntrlArr = [3, 3];
% выбор формы ДН для алгоритма управления шириной ДН (antPattCntrl = 3)
% 0 - окно Гаусса
% 1 - окно приподнятого косинуса
% 2 - прямоугольное окно
win_typeArr = [2, 0];
antTypeCmt = ["URA", "ULA", "UCA"];
Zaa = [];
for aa=1:length(antPattCntrlArr)
    Zst = [];
    antPattCntrl = antPattCntrlArr(aa);
    win_type = win_typeArr(aa);
    switch win_type
        case 0
            sclArr = 0.1:0.5:5.5; % окно Гаусса 
        case 1
            sclArr = 0:0.1:0.9;   % окно приподн. косинуса; не моделируется
        case 2
            sclArr = 0.5:0.5:5.5; % прямоугольное окно
    end
    for s=stdCoordsArr
        Zn = [];
        stdCoords = s;
        for nel=NellArr
            Zsc = [];
            Nel = nel;
            for scl=sclArr
                rng('default');
                fprintf('win %i std=%.1f, Nel=%i, Scl=%.2f\n', ...
                    win_type, s, nel, scl);
                [X,Y,Z] = lab_beam_shape_fcn(antType,Nel,stdCoords,...
                            antPattCntrl,win_type,scl);
                Zsc = [Zsc, mean(Z(:))];
            end
            Zn = [Zn;Zsc];
        end
        Zst{end+1} = Zn;
    end
    Zaa{end+1} = Zst;
end

% построение графиков
for aa=1:length(antPattCntrlArr)
    figure;
    ZZmS = Zaa{aa};
    [X,Y] = meshgrid(NellArr, sclArr);
    for kk=1:size(ZZmS,2)
        subplot(2,2,kk)
        surf(X, Y, ZZmS{kk}.', 'FaceColor', 'interp', 'EdgeColor','none');
        c1 = colorbar; c1.Label.String = 'SIR_{avg}, дБ';
    grid on; xlabel('N'); ylabel('s'); view([0, 90]); axis tight
    legend(sprintf('RMSE = %d м',stdCoordsArr(kk)),'Location','southeast');
    title(sprintf('RMSE = %d м',stdCoordsArr(kk)));
    end
    if win_typeArr(aa) == 0
        sgtitle(sprintf('ДН Гаусса; %s', antTypeCmt(antType)));
    elseif win_typeArr(aa) == 2
        sgtitle(sprintf('Прямоугольная ДН; %s', antTypeCmt(antType)));
    end
end