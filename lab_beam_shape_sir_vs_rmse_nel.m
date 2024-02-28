% построение графиков зависимости SIR и & delta SIR от числа элементов АР и
% СКО оценки координат для сценариев управления: 1) формой луча с 
% максимумом и нулем ДН; 2) шириной луча с ДН Гаусса; 3) шириной луча с 
% прямоугольной ДН Гаусса; зависимости delta SIR построены относительно
% зависимости SIR при управлении только максимумом ДН; зависимости SIR для
% управления только максимумом ДН построены пунктиром точек
clear all; close all; clc;
stdCoordsArr = [1, 5, 10, 15];  % СКО оценки координат UE в метрах
NellArr = 2:2:30;               % число АЭ в одном измерении
% выбор типа антенной решетки
% 1 - planar or uniform rectangural antenna array (URA), планарная АР
% 2 - uniform linear antenna array (ULA), линейная АР
% 3 - uniform circular antenna array (UCA); круговая АР: antPattCntrl=0,1,2
antType = 1;
% вектор моделей управления формой и шириной луча при ДО:
%   0 - управление формой луча: максимумом ДН 
%   1 - управление формой луча: максимумом и нулем ДН 
%   2 - адаптивное управление формой/шириной ДНА (не зависит от RMSE)
%   3 - управление шириной луча 
antPattCntrlArr = [0, 1, 3, 3];
antPattCntrlCmt = ["управление максимумом ДН", ...
                   "управление максимумом и нулем ДН", ...
                   "управление шириной прямоугольной ДН",...
                   "управление шириной ДН Гаусса"];
antTypeCmt = ["URA", "ULA", "UCA"];
% выбор формы ДН для алгоритма управления шириной луча (antPattCntrl = 3)
% 0 - окно Гаусса
% 1 - окно приподнятого косинуса
% 2 - прямоугольное окно
win_typeArr = [0, 0, 2, 0];

Zaa = [];
for aa=1:length(antPattCntrlArr)
    Zst = [];
    antPattCntrl = antPattCntrlArr(aa);
    win_type = win_typeArr(aa);
    for s=stdCoordsArr
        Zn = [];
        stdCoords = s;
        for nel=NellArr
            Zsc = [];
            Nel = nel;
            switch win_type
                case 0
                    sc = 3.79*exp(-0.126*stdCoords);
                case 1
                    sc = 0.1;
                case 2
                    sc = -0.2*stdCoords + 4.65;
            end
            rng('default');
            fprintf('alg %i rmse=%5.1f, nel=%3i, scl=%4.2f\n',...
                antPattCntrl, s, nel, sc);
            [X,Y,Z] = lab_beam_shape_fcn(antType,Nel,stdCoords,...
                                            antPattCntrl,win_type);
            Zsc = [Zsc, mean(Z(:))];
            Zn = [Zn;Zsc];
        end
        Zst{end+1} = Zn;
    end
    Zaa{end+1} = Zst;
end
% построение графиков
ZZmS_1 = Zaa{1};
for aa=2:length(antPattCntrlArr)
    ZZmS = Zaa{aa};
    figure; hold on; grid on; axis tight;
    for kk=1:size(ZZmS,2)
        plot(NellArr, ZZmS_1{kk}, '-', 'linewidth',2); hold on;
    end
    set(gca,'ColorOrderIndex',1);
    for kk=1:size(ZZmS,2)
        plot(NellArr, ZZmS{kk}, '--', 'linewidth',2); hold on;
    end
    xlabel('N (в одном измерении для URA)'); ylabel('SIR_{avg}, дБ');
    legFrst = ['RMSE = ',num2str(stdCoordsArr(1)), ' м'];
    legend([legFrst, strcat(string(stdCoordsArr(2:end)),' м')],...
        'NumColumns',3, 'Location', 'northwest');
    title([antPattCntrlCmt(aa)]);
    
    figure; hold on; grid on; axis tight;
    for kk=1:size(ZZmS,2)
        plot(NellArr, ZZmS{kk}-ZZmS_1{kk}, '-.', 'linewidth',2); hold on;
    end
    xlabel('N (в одном измерении для URA)'); ylabel('\DeltaSIR_{avg}, дБ');
    legFrst = ['RMSE = ',num2str(stdCoordsArr(1)), ' м'];
    legend([legFrst, strcat(string(stdCoordsArr(2:end)),' м')],...
        'NumColumns',3, 'Location', 'northwest');
    title([antPattCntrlCmt(aa)]);
end