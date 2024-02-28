clear all; close all; clc;
c = physconst('LightSpeed');
anim = 1;       % 1 - вкл. анимацию работы
stdCoords = 0;  % СКО оценки координат UE по осям x, y, z в метрах
f = 30e9;       % несущая в диапазоне ММВ, Гц
lamb = c/f;     % длина волны, м
da = 0.5*c/f;   % расстояние между элементами АР
snrThr = 10;    % пороговое отношения сигнал/помеха для отображения карты
useAntUE = 0;   % использовать ДН АР на UE (только режим antPattCntrl=0)
% выбор алгоритма управления ДН на NB:
% 0 - управление максимумом ДН
% 1 - управление максимумом и нулем ДН
% 2 - адаптивное управление ДН
% 3 - управление шириной ДН
antPattCntrl = 0;
% выбор формы ДН для алгоритма управления шириной луча (antPattCntrl = 3)
% 0 - окно Гаусса
% 1 - окно приподнятого косинуса
% 2 - прямоугольное окно
win_type = 2;
backLobe = 1; % подавление обратного лепестка ДН (логично только для URA)
% отображение ДН в 2D при выставленном флаге анимации (anim = 1)
antPlot2D = 1; % (рекомендуется выставлять в 1)
% Шаг отрисовки ДН при anim = 1
if (antPlot2D == 0)
    angStep = 5;
else
    angStep = 1;
end
% период процедуры диаграммообразования, с;
Ta = 0; % при Ta = 0 ДН формируется на каждом шаге расчета
% выбор типа антенной решетки
% 1 - планарная АР
% 2 - линейная АР
% 3 - круговая АР (поддержаны только режимы ДН antPattCntrl=0,1,2)
antType = 2;
Nel = 10;                               % число АЭ в одном измерении
antElPos = createAnt(antType, Nel, da); % формирование АР
NelFull = size(antElPos, 1);            % полное число АЭ
% выбор сценария:
% 1 - gNB расположены на одной стороне относительно траектории движения UE
% 2 - gNB расположены по разные стороны относительно траектории движения UE
[gNB, ueNode, d, T, v] = createScenarion(1, 10);
Nd = length(d);

Nnb = length(gNB);    % число gNB
Nue = length(ueNode); % число UE
% число точек расчета (число точек координат траектории UE)
N = length(ueNode(1).Trajectory(:,1));
trajArray = [ueNode.Trajectory];       % массив координат UE [N x 3*Nue]
gNBcoords = [gNB(:).Coords].';         % массив координат gNB [Nnb x 3]

if (anim == 1)
    antPattScl = 30.0;  % коэфф. масштабирования для отображения ДН gNB
    antPattSclUe = 5.0; % коэфф. масштабирования для отображения ДН UE
    fg = figure(5); fg.WindowState = 'maximized'; grid on; hold on;
    % отображение ДН gNB (параметры см. в antPattPlot)
    gNBptrnPlot = gobjects(1, 2);
    for i=1:Nnb
        [x, y, z] = antPattPlot(antElPos, f, gNB(i), ...
                                angStep, antPattScl, backLobe, antPlot2D);
        if (antPlot2D == 0)
            gNBptrnPlot(i) = surf(x+gNB(i).Coords(1),...
                            y+gNB(i).Coords(2),...
                            z+gNB(i).Coords(3), 'FaceColor', '#4DBEEE');
        else
            gNBptrnPlot(i) = plot(x+gNB(i).Coords(1),...
                            y+gNB(i).Coords(2),'Color', '#4DBEEE');
        end
    end
    
    % отображение ДН UE (параметры см. в antPattPlot)
    if (useAntUE == 1)
        ueptrnPlot = gobjects(1, 2);
        for i=1:Nue
            [x, y, z] = antPattPlot(antElPos, f, ueNode(i), ...
                angStep, antPattSclUe, backLobe);
            ueptrnPlot(i) = surf(x+ueNode(i).Trajectory(1,1),...
                            y+ueNode(i).Trajectory(1,2),...
                            z+ueNode(i).Trajectory(1,3), ...
                            'FaceColor', '#4DBEEE');
        end
    end

    % отображение положения UE
    uePlot = gobjects(1, 2);
    ueText = gobjects(1, 2);
    for i=1:Nue
        uePlot(i) = plot3(ueNode(i).Trajectory(1,1),...    
                          ueNode(i).Trajectory(1,2),...
                          ueNode(i).Trajectory(1,3), '^', ...
                          'MarkerSize', 10);
        ueText(i) = text(ueNode(i).Trajectory(1,1), ...
                         ueNode(i).Trajectory(1,2), ...
                         ueNode(i).Trajectory(1,3), ...
                         sprintf('UE_{%i}', i), 'FontSize', 14, ...
                         'Color', '#A2142F');
    end
    ueDirPlot  = gobjects(1, 2);
    ueDirPlot2 = gobjects(1, 2);
    indSnoi = [2,1];
    for i=1:Nnb
        % отображение вектора направления от gNB на UE
        ueDirPlot(i)=plot3([gNB(i).Coords(1);ueNode(i).Trajectory(1,1)],...
                           [gNB(i).Coords(2);ueNode(i).Trajectory(1,2)],...
                           [gNB(i).Coords(3);ueNode(i).Trajectory(1,3)],...
                           'Color', '#76AB2F');
        % отображение вектора направления от gNB на соседней UE
        ueDirPlot2(i) = plot3(...
            [gNB(i).Coords(1);ueNode(indSnoi(i)).Trajectory(1,1)],...
            [gNB(i).Coords(2);ueNode(indSnoi(i)).Trajectory(1,2)],...
            [gNB(i).Coords(3);ueNode(indSnoi(i)).Trajectory(1,3)], ...
            'Color', '#D95319');
    end
    ueDirPlot(2).LineStyle = '--';
    ueDirPlot2(2).LineStyle = '--';
    for i=1:Nnb
        text(gNB(i).Coords(1), gNB(i).Coords(2), gNB(i).Coords(3)+5, ...
            sprintf('gNB_{%i}', i), 'FontSize', 16, 'Color', '#A2142F');
    end
    xlabel('x, m'); ylabel('y, m'); axis equal;
    axis([-5, 155, -5, 65, 0, 30]); view([0, 90]);
end

% инициализация массивов для хранения углов отправки от каждой gNB 
% к каждой UE (нужно для устранении ошибки чтения несуществующего 
% массива при некоторых настройках модели)
azAng = zeros(Nue, Nnb);
elAng = zeros(Nue, Nnb);
azAngUE = zeros(Nue, Nnb);
elAngUE = zeros(Nue, Nnb);

%% ЦИКЛ ПО ЧИСЛУ ТОЧЕК РАСЧЕТА
for i=1:76 % цикл по числу точек расчета
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
            % вектор, задающий направление из gNB в UE в глобальной СКx,y,z 
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
                azAngT(j, n) = rad2deg(atan2(dirVectT(2), dirVectT(1)));
                elAngT(j, n) = rad2deg(atan2(dirVectT(3), ...
                    sqrt(sum(dirVectT(1:2).^2))));
                % расчет углов ухода для АР UE
                if ( useAntUE == 1)
                    % вектор, задающий направление из UE в gNB в локальной
                    % системе координат АР UE,т.е. с учетом положения АР UE
                    dirVect = -ueNode(j).AntOrient.'*diffCoord(n,:).';
                    % расчет углов ухода от j-й UE к n-й gNB
                    azAngUE(j, n) = rad2deg(atan2(dirVect(2), dirVect(1)));
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
                cVect(:,j) = ...
                 getAntPatternSteer(antElPos, f, azAng(j, n), elAng(j, n));
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
                    sqrt(sum((ueCoordsiErr(j,:).' - gNB(n).Coords).^2))];
                end
                % расчет вектора направляющих коэффициентов
                % АР j-й gNB работающей с n-й gNB
                if (n == ueNode(j).servgNB && useAntUE == 1)
                    ueNode(j).Steer = getAntPatternSteer(antElPos, f, ...
                        azAngUE(j,n), elAngUE(j,n));
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
                            scl = 0.45;
                        case 1
                            scl = 0.1;
                        case 2
                            scl = 2.2;
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
                gUE = getAntPatternG(antElPos, f, ...
                    azAngUE(n), elAngUE(n), ueNode(1).Steer, backLobe).^2;
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
    % обновление отображения положения UE и ДН АР gNB/UE 
    % для текущего шага расчета
    if (anim == 1)
        for ip=1:Nnb
            % обновление ДН gNB
            [x, y, z] = antPattPlot(antElPos, f, ...
                gNB(ip), angStep, antPattScl, backLobe, antPlot2D);
            gNBptrnPlot(ip).XData = x + gNB(ip).Coords(1);
            gNBptrnPlot(ip).YData = y + gNB(ip).Coords(2);
            if (antPlot2D == 0)
                gNBptrnPlot(ip).ZData = z + gNB(ip).Coords(3);
            end
            % обновление ДН UE
            if (useAntUE == 1)
                [x, y, z] = antPattPlot(antElPos, f, ...
                    ueNode(ip),angStep, antPattSclUe, backLobe, antPlot2D);
                ueptrnPlot(ip).XData = x + ueCoordsi(ip,1);
                ueptrnPlot(ip).YData = y + ueCoordsi(ip,2);
                if (antPlot2D == 0)
                    ueptrnPlot(ip).ZData = z + ueCoordsi(ip,3);
                end
            end
            % обновление вектора направления от gNB на UE
            ueDirPlot(ip).XData = [gNB(ip).Coords(1); ueCoordsi(ip,1)];
            ueDirPlot(ip).YData = [gNB(ip).Coords(2); ueCoordsi(ip,2)];
            ueDirPlot(ip).ZData = [gNB(ip).Coords(3); ueCoordsi(ip,3)];
            % обновление вектора направления от gNB на  соседнюю UE
            ueDirPlot2(ip).XData = [gNB(ip).Coords(1); ...
                ueCoordsi(indSnoi(ip),1)];
            ueDirPlot2(ip).YData = [gNB(ip).Coords(2); ...
                ueCoordsi(indSnoi(ip),2)];
            ueDirPlot2(ip).ZData = [gNB(ip).Coords(3); ...
                ueCoordsi(indSnoi(ip),3)];
            % обновление положения UE
            uePlot(ip).XData = ueCoordsi(ip,1);
            uePlot(ip).YData = ueCoordsi(ip,2);
            uePlot(ip).ZData = ueCoordsi(ip,3);
            ueText(ip).Position = [ueCoordsi(ip,1)+2, ...
                                   ueCoordsi(ip,2), ...
                                   ueCoordsi(ip,3)+5];
        end  % for ip=1:Nnb    
        pause(0.01)
    end % if (anim == 1)
end % for i=1:N

%%
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
    figure(1); surf(X, Y, Z, 'FaceColor', 'interp', 'EdgeColor','none');
    grid on; xlabel('x, м'); ylabel('y, m'); view([0, 90]);
    c1 = colorbar; c1.Label.String = 'SIR, дБ'; 
    % карта точек положения 1-й UE, в которых отношение сигнал/помеха 
    % превышает заданный порог snrThr
    figure(2); surf(X, Y, double(Z>snrThr), 'EdgeColor','none');
    grid on; xlabel('x, м'); ylabel('y, м'); view([0, 90]);
    colormap(winter(2)); c3 = colorbar;
    c3.Label.String = sprintf('SIR > %.0f дБ', snrThr);
    c3.Ticks = [0, 1]; view([0, 90]);   
    % карта отношения сигнал/помеха с отображением положения 
    % и ориентации АР gNB
    figure(3); surf(X, Y, Z, 'FaceColor', 'interp', 'EdgeColor','none');
    grid on; hold on; 
    for n=1:Nnb
        absAntCoord = (gNB(n).AntOrient*[0, -1, 0; 0, 1, 0].'*3 + ...
            gNB(n).Coords).';
        plot3(absAntCoord(:,1), absAntCoord(:,2), absAntCoord(:,3), ...
            'Color', '#ECB01F', 'LineWidth', 3);
        text(gNB(n).Coords(1), gNB(n).Coords(2)*1.06, gNB(n).Coords(3), ...
            sprintf('Ant gNB %i', n));
    end
    xlabel('x, м'); ylabel('Y, м'); axis equal; view([0, 90]);
    c2 = colorbar; c2.Label.String = 'SIR, дБ';
end