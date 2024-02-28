clear all; close all; clc;
stdCoords = 0;  % СКО оценки координат UE по осям x, y, z в метрах
Nel = 10;       % число АЭ в одном измерении
c = physconst('LightSpeed');
f = 30e9;     % несущая в диапазоне ММВ, Гц
lamb = c/f;   % длина волны, м
da = 0.5*c/f; % расстояние между элементами АР
snrThr = 10;  % пороговое отношения сигнал/помеха для отображения карты, дБ
useAntUE = 0; % использовать ДН АР на UE (только режим antPattCntrl=0)
backLobe = 1; % подавление обратного лепестка ДН (логично только для URA)
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
antPattCntrlArr = [0, 1, 2, 3, 3];
% выбор формы ДН для алгоритма управления шириной луча (antPattCntrl = 3)
% 0 - окно Гаусса
% 1 - окно приподнятого косинуса
% 2 - прямоугольное окно
win_typeArr = [0, 0, 0, 0, 2];
antPattCntrlCmt = ["управление максимумом ДН", ...
                   "управление максимумом и нулем ДН", ...
                   "адаптивное управление ДН", ...
                   "управление шириной ДН Гаусса", ...
                   "управление шириной прямоугольной ДН"];
antTypeCmt = ["URA", "ULA", "UCA"];
figNumber = 0;
for aa=1:length(antPattCntrlArr)
    antPattCntrl = antPattCntrlArr(aa); % алгоритм управления лучом
    win_type = win_typeArr(aa); % форма ДН для управления шириной луча
    % период процедуры диаграммообразования, с;
    Ta = 0; % при Ta = 0 ДН формируется на каждом шаге расчета
    antElPos = createAnt(antType, Nel, da); % формирование АР
    NelFull = size(antElPos, 1);            % полное число АЭ
    % выбор сценария:
    % 1 - gNB расположены на одной стороне относительно траектории UE
    % 2 - gNB расположены по разные стороны относительно траектории UE
    [gNB, ueNode, d, T, v] = createScenarion(1, 0:0.1:10);
    Nd = length(d);
    Nnb = length(gNB);    % число gNB
    Nue = length(ueNode); % число UE
    % число точек расчета (число точек координат траектории UE)
    N = length(ueNode(1).Trajectory(:,1));
    trajArray = [ueNode.Trajectory];      % массив координат UE [N x 3*Nue]
    gNBcoords = [gNB(:).Coords].';        % массив координат gNB [Nnb x 3]
    % инициализация массивов для хранения углов отправки от каждой gNB 
    % к каждой UE (нужно для устранении ошибки чтения несуществующего 
    % массива при некоторых настройках модели)
    azAng = zeros(Nue, Nnb);
    elAng = zeros(Nue, Nnb);
    azAngUE = zeros(Nue, Nnb);
    elAngUE = zeros(Nue, Nnb);
    % ЦИКЛ ПО ЧИСЛУ ТОЧЕК РАСЧЕТА
    for i=1:N % цикл по числу точек расчета
        % массив координат всех UE для i-й точки расчета
        ueCoordsi = reshape(trajArray(i, :).', 3, Nue).';
        % внесение ошибки в оценку координат UE согласно stdCoords
        ueCoordsiErr = ueCoordsi;
        ueCoordsiErr(:,1:2) = ueCoordsiErr(:,1:2) + stdCoords*randn(2,2);
        if (mod(i, round(Ta/T)) == 1 || Ta == 0)
            % инициализация массивов для хранения углов ухода 
            % от каждой gNB к каждой UE
            azAng = zeros(Nue, Nnb);  % угол ухода по азимуту (AOD)
            elAng = zeros(Nue, Nnb);  % угол ухода по углу места (AOD)
            azAngT = zeros(Nue, Nnb); % истинный угол ухода по азимуту (AOD)
            elAngT = zeros(Nue, Nnb); % истинный угол ухода по углу места (AOD)
            % инициализация массивов для хранения углов прихода 
            % для каждой UE от каждой gNB (используется при наличии АР на UE)
            azAngUE = zeros(Nue, Nnb);
            elAngUE = zeros(Nue, Nnb);
            for j=1:Nue % цикл по числу UE
                % вектор, задающий направление из gNB в UE в глобальной СК x,y,z 
                diffCoord = ueCoordsiErr(j,:) - gNBcoords;
                diffCoordT = ueCoordsi(j,:) - gNBcoords;
                for n=1:Nnb % цикл по числу gNB
                    % вектор, задающий направление из gNB в UE в локальной 
                    % системе координат АР gNB (т.е. с учетом положения АР gNB)
                    dirVect = gNB(n).AntOrient.'*diffCoord(n,:).';
                    dirVectT = gNB(n).AntOrient.'*diffCoordT(n,:).';
                    % расчет углов ухода от n-й gNB к j-й UE
                    azAng(j, n) = rad2deg(atan2(dirVect(2), dirVect(1)));
                    elAng(j, n) = rad2deg(atan2(dirVect(3), ...
                    sqrt(sum(dirVect(1:2).^2))));
                    azAngT(j, n) = rad2deg(atan2(dirVectT(2),dirVectT(1)));
                    elAngT(j, n) = rad2deg(atan2(dirVectT(3), ...
                    sqrt(sum(dirVectT(1:2).^2))));
                    % расчет углов ухода для АР UE
                    if ( useAntUE == 1)
                        % вектор, задающий направление из UE в gNB в локальной
                        % системе координат АР UE,т.е. с учетом положения АР UE
                        dirVect = -ueNode(j).AntOrient.'*diffCoord(n,:).';
                        % расчет углов ухода от j-й UE к n-й gNB
                        azAngUE(j, n) = rad2deg(atan2(dirVect(2), ...
                                dirVect(1)));
                        elAngUE(j, n) = rad2deg(atan2(dirVect(3), ...
                                sqrt(sum(dirVect(1:2).^2))));
                    end
                end
            end
            % для адаптивной схемы, т.к. она не использует информацию 
            % о местоположении UE
            if (antPattCntrl == 2)
                azAng = azAngT;
                elAng = elAngT;
            end
            for n=1:Nnb % цикл по числу gNB
                cVect = zeros(NelFull, Nue);
                gVect = zeros(1, Nue);
                for j=1:Nue % цикл по числу UE      
                    % расчет направляющего вектора 
                    % фазового распределения АР n-й gNB к j-й UE
                    cVect(:,j) = getAntPatternSteer(antElPos, f,...
                                    azAng(j, n), elAng(j, n));
                    % сохранения данных о обслуживаемой UE (номер обслуживающей 
                    % gNB указывается в параметре servgNB для каждой UE)
                    if (n == ueNode(j).servgNB)
                        % расчет вектора направляющих коэффициентов
                        % АР n-й gNB обслуживающей j-ю UE  
                        gNB(n).Steer = getAntPatternSteer(antElPos, f, ...
                            azAng(j, n), elAng(j, n))/NelFull;
                        % запись 1 в векторе g (см. формулу (10)) по индексу
                        % j-й UE для выставления максимума в этом направлении
                        gVect(j) = 1;
                        % сохранение информации о местоположении UE в полярных
                        % координатах (азимут, угол места, расстояние)
                        gNB(n).UEPolarCoord = [azAng(j, n), elAng(j, n),...
                        sqrt(sum((ueCoordsiErr(j,:).'-gNB(n).Coords).^2))];
                    end
                    % расчет вектора направляющих коэффициентов
                    % АР j-й gNB работающей с n-й gNB
                    if (n == ueNode(j).servgNB && useAntUE == 1)
                        ueNode(j).Steer = getAntPatternSteer(antElPos,...
                                f, azAngUE(j,n), elAngUE(j,n));
                    end
                end
                % расчет вектора направляющих коэфф. АР n-й gNB обслуживающей
                % j-ю UE при схеме управления ДН отличной от antPattCntrl=0
                switch antPattCntrl
                    % LCMV
                    case 1
                        % расчет вектора коэффициентов по алгоритму LCMV для
                        % управления максимумом и нулем ДН
                        gNB(n).Steer = cVect*pinv(cVect'*cVect)*gVect';
                        % нормировка коэффициентов
                        gNB(n).Steer = ...
                            gNB(n).Steer/max(abs(gNB(n).Steer))/NelFull;
                    % SMI
                    case 2
                        K = 100; % длина преамбулы
                        % направляющий вектор АР для полезного сигнала SOI
                        sv_s = cVect(:,[ueNode.servgNB] == n);
                        % направляющий вектор АР для помехи SNOI
                        sv_j = cVect(:,~([ueNode.servgNB] == n));
                        % симуляция приема сигнала и помехи с заданных
                        % направлений
                        rx_sum = ones(K,1)*sv_s.' + ...
                            (randn(K,1) + 1i*randn(K,1))/sqrt(2)*sv_j.';
                        % корреляционная матрица сигналов АЭ
                        Rxx = rx_sum'*rx_sum;
                        % корреляционный вектор сигналов АЭ 
                        % с опорной преамбулой
                        rxd = rx_sum'*ones(K,1);
                        % расчет вектора коэфф. по алгоритму SMI 
                        % для подавления помехи
                        gNB(n).Steer = conj(pinv(Rxx)*rxd);
                    % beamshaping
                    case 3
                        % множитель ширины ДН
                        switch win_type
                            case 0
                                scl = 3.79*exp(-0.126*stdCoords);
                            case 1
                                scl = 0.1;
                            case 2
                                scl = -0.2*stdCoords + 4.65;
                        end
                        % ширина ДН, рассчитанная по СКО оценки координат UE
                        BW = 2*atan2d(stdCoords,gNB(n).UEPolarCoord(3));
                        % расчет вектора весовых коэфф. АР
                        w = beamshapingWeight(win_type, BW, ...
                            gNB(n).UEPolarCoord(1), Nel, scl);
                        % применять коэфф. w если требуемая ширина ДН:
                        % не меньше, чем минимальная теоретическая ширина
                        % равная 0.891*lamb/Nel/da; амплитуда суммы коэфф.>0.1
                        if ~(BW*scl <= rad2deg(0.891*lamb/Nel/da) ||...
                                any(isnan(w)) || sum(w) < 0.1)
                            if (antType == 1)
                                % расчет вектора весовых коэффициентов
                                % планарной АР для вертикальных АЭ
                                antElPosElev = [zeros(Nel,1), ...
                                    antElPos(1:Nel,2), zeros(Nel,1)];
                                wEl = getAntPatternSteer(antElPosElev, ...
                                    f, gNB(n).UEPolarCoord(2), 0)/Nel;
                                % вектор весовых коэффициентов планарной АР,
                                % полученный перемножением горизонтальных и 
                                % вертикальных коэффициентов
                                w = w*wEl.'/Nel;
                                w = w(:);
                            end                        
                            gNB(n).Steer = w;
                        end
                end % switch antPattCntrl
            end % for n=1:Nnb     
        end % if (mod(i, round(Ta/T)) == 1 || Ta == 0)
        % массив для временного хранения значений 
        % принимаемой мощности на UE от каждой gNB
        gNBpwr = zeros(Nnb, 1);
        % расчет отношения принимаемой мощности от обслуживающей gNB 
        % к мощности, принимаемой от соседней gNB для каждой UE
        for j=1:Nue % цикл по числу UE
            for n=1:Nnb % цикл по числу gNB
                % расчет коэфф усиления приемной антенны UE
                if ( useAntUE == 1 && j == 1)
                    % КУ при наличии АР на UE, луч которой 
                    % направлен на обслуживающую gNB
                    gUE = getAntPatternG(antElPos, f, azAngUE(n),...
                        elAngUE(n), ueNode(1).Steer, backLobe).^2;
                else
                    % КУ без АР на UE
                    gUE = 1;
                end
                % расчет мощности принимаемой на j-й UE от n-й gNB с учетом
                % диаграммообразования на UE и gNB (без учета дальности)
                gNBpwr(n) = gUE*getAntPatternG(antElPos, f, ...
                    azAngT(j, n), elAngT(j, n), gNB(n).Steer, backLobe).^2;
            end
            % расчет расстояния от j-й UE до каждой gNB
            diffCoord = ueCoordsi(j,:) - gNBcoords;
            distSpace = sqrt(sum(diffCoord.^2,2));
            % расчет мощности принимаемой на j-й UE от каждой gNB с учетом
            % дальности (потери рассчитываются по модели затухания в свободном
            % пространстве)
            gNBpwr = pow2db(gNBpwr) - fspl(distSpace,c/f);
            % расчет отношения принимаемой мощности от обслуживающей gNB 
            % к мощности принимаемой от соседней gNB для j-й UE
            ueNode(j).SNR(i) = gNBpwr(ueNode(j).servgNB) - ...
                    sum(gNBpwr(1:end ~= ueNode(j).servgNB));
        end
    end % for i=1:N 
    % подготовка массивов координат и массива значений 
    % отношения сигнал/помеха 1-й UE для отображения
    X = reshape(ueNode(1).Trajectory(:,1), [], Nd).';
    Y = reshape(ueNode(1).Trajectory(:,2), [], Nd).';
    Z = reshape(ueNode(1).SNR, [], Nd).';
    if size(Z, 1) == 1
        figure(1); plot(X, Z); grid on;
        xlabel('x, м'); ylabel('SIR, дБ');
    else
        % карта отношения сигнал/помеха в каждой точке положения 1-й UE
        figure; surf(X, Y, Z, 'FaceColor', 'interp', 'EdgeColor','none');
        grid on; xlabel('X, м'); ylabel('Y, м'); view([0, 90]);
        c1 = colorbar; c1.Label.String = 'SIR, дБ'; 
        title([antPattCntrlCmt(aa), antTypeCmt(antType)]);
        % карта точек положения 1-й UE, в которых отношение сигнал/помеха 
        % превышает заданный порог snrThr
        figure; surf(X, Y, double(Z>snrThr), 'EdgeColor','none');
        grid on; xlabel('x, м'); ylabel('y, м'); view([0, 90]);
        colormap(winter(2)); c3 = colorbar;
        c3.Label.String = sprintf('SIR > %.0f дБ', snrThr);
        c3.Ticks = [0, 1]; view([0, 90]);   
        title([antPattCntrlCmt(aa), antTypeCmt(antType)]);
    end
    figNumber = figNumber + 1;
end % for aa=1:length(antPattCntrlArr)