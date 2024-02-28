clear all; close all; clc;
c = physconst('LightSpeed');
% выбор типа антенной решетки:
% 1 - planar or uniform rectangural antenna array (URA), планарная АР
% 2 - uniform linear antenna array (ULA), линейная АР
% 3 - uniform circular antenna array (UCA); круговая АР: не поддерживается
antType = 1;
Nel = 20;       % число АЭ в одном измерении
f = 30e9;       % несущая в диапазоне ММВ, Гц
lamb = c/f;     % длина волны, м
da = 0.5*c/f;   % расстояние между элементами АР
antElPos = createAnt(antType, Nel, da); % формирование АР
NelFull = size(antElPos, 1);            % общее число АЭ в АР
% выбор формы ДН для алгоритма управления шириной ДН (antPattCntrl = 3)
% 0 - окно Гаусса
% 1 - окно приподнятого косинуса
% 2 - прямоугольное окно
win_typeArr = [2, 0];
antTypeCmt = ["URA", "ULA", "UCA"];
ueRxPwrPlt = [];
figNumber = 1;
stdCoordsArr=[10];
for s=stdCoordsArr
    subFigNumber = 1;
    for ww=win_typeArr
        stdCoords = s; % СКО оценки координат UE, м
        N = 100;       % число точек расчета
        HgNB = 0;      % высота подвеса АР gNB
        Due = 50;      % расстояние между gNB и UE на плоскости
        % структура параметров gNB
        gNB = createNB([0, 0, HgNB], [0, 0]);
        gNB.Steer = zeros(NelFull, 2);
        ueRxPwr = zeros(2, N);
        ueCoord = [Due, 0, 0];         % координаты UE
        gNBcoords = [gNB(:).Coords].'; % массив координат gNB
        % gNB-UE 3D расстояние
        distSpaceT = sqrt(sum((gNB.Coords-ueCoord.').^2)); 
        BW = 2*atan2d(stdCoords,distSpaceT);        % ширина ДН
        stAng = 0;                                  % направление максимума
        % вектор, задающий направление из gNB в UE в глобальной СК x,y,z
        diffCoord = ueCoord - gNBcoords;
        % вектор, задающий направление из gNB в UE 
        % в локальной СК АР gNB (т.е. с учетом положения АР gNB)
        dirVect = gNB.AntOrient.'*diffCoord.';
        % расчет углов ухода AOD от gNB к UE 
        azAng = rad2deg(atan2(dirVect(2), dirVect(1)));
        elAng = rad2deg(atan2(dirVect(3), sqrt(sum(dirVect(1:2).^2))));
        % расчет вектора направляющих коэффициентов АР gNB
        gNB.Steer(:,1)=getAntPatternSteer(antElPos,f,azAng,elAng)/NelFull;
        scl = 1.5; % масштабирующий множитель ширины ДН
        % расчет вектора весовых коэфф. АР
        [w, azAngP, antPattP] = beamshapingWeight(ww, BW, stAng, Nel, scl);
        % применять коэффициенты w если требуемая ширина ДН не меньше, 
        % чем минимальная теоретическая ширина 0.891*lamb/Nel/da
        if (BW*scl <= rad2deg(0.891*lamb/Nel/da) || any(isnan(w)))
            w = gNB.Steer(:,1);
        else
            if (antType == 1)
                % расчет вектора весовых коэфф. планарной АР 
                % для вертикальных АЭ
                w = repmat(w, Nel, 1)/Nel;
                w = w(:);
            end
        end
        gNB.Steer(:,2) = w;
        rng('default');
        for i=1:N % цикл по числу точек расчета
            % внесение ошибки в оценку координат UE по stdCoords
            ueCoordErr = ueCoord;
            ueCoordErr(1:2) = ueCoordErr(1:2) + ...
                stdCoords*randn(size(ueCoord(1:2)));
            % вектор, задающий направление из gNB в UE в глобальной 
            % системе координат x,y,z с учетом ошибки координат UE
            diffCoordT = ueCoordErr - gNBcoords;
            % вектор, задающий направление из gNB в UE в системе
            % координат АР gNB (т.е. с учетом положения АР gNB)
            dirVectT = gNB.AntOrient.'*diffCoordT.';
            % расчет углов отправки от gNB к UE
            azAngT = rad2deg(atan2(dirVectT(2), dirVectT(1)));
            elAngT=rad2deg(atan2(dirVectT(3),sqrt(sum(dirVectT(1:2).^2))));
            % расчет принимаемой мощности от обслуживающей gNB 
            % с учетом диаграммообразования на gNB (без учета 
            % дальности) без/с управления шириной ДН  
            gNBpwr = [getAntPatternG(antElPos, f, ...
                        azAngT, elAngT, gNB.Steer(:,1), 0).^2;...
                        getAntPatternG(antElPos, f, ...
                        azAngT, elAngT, gNB.Steer(:,2), 0).^2];
            % расчет расстояния от UE до gNB
            diffCoord = ueCoordErr - gNBcoords;
            distSpace = sqrt(sum(diffCoord.^2,2));
            % расчет мощности принимаемой UE от gNB с учетом дальности; 
            % потери рассчитываются по модели затухания в свободном
            % пространстве для случая без/с управлением шириной ДН
            gNBpwr(isnan(gNBpwr)) = gNBpwr(1);
            gNBpwr = pow2db(gNBpwr) - fspl(distSpace,c/f);
            ueRxPwr(:, i) = gNBpwr;
        end % for i=1:N 
        % расчет ДН для случая без/с управление шириной луча HPBW
        azAngPatt = -90:0.1:89;
        g = zeros(1, length(azAngPatt));
        gDef = zeros(1, length(azAngPatt));
        for i=1:length(azAngPatt)
            g(i) = getAntPatternG(antElPos, f, azAngPatt(i), 0, w, 0);
            gDef(i) = getAntPatternG(antElPos, f, ...
            azAngPatt(i), 0, ones(size(w))/NelFull, 0);
        end
        gNorm = g/max(g);
        % расчет ширины луча HPBW по уровню - 3 дБ
        [~,ind] = min(abs(gNorm - 1/sqrt(2)));
        hpbw = 2*abs(azAngPatt(ind));
        % расчет координат кривых ДН для отображения
        % повернутого на 90 градусов для лучшей визуализации
        alphPatt = azAngPatt + 90;
        xPatt = cosd(alphPatt).*gNorm*(distSpaceT);
        yPatt = sind(alphPatt).*gNorm*(distSpaceT);
        xPattDef = cosd(alphPatt).*gDef*(distSpaceT);
        yPattDef = sind(alphPatt).*gDef*(distSpaceT);        
        % отображение ДН на карте вероятности положения UE
        figure(figNumber); 
        subplot(1,2,subFigNumber); 
        prbA = [0.1, 0.5, 0.9];
        [~, ~, ppR] = get_prob(ueCoord(2), ueCoord(1), stdCoords, prbA);
        
        ylimnin=-2*stdCoords; ylimnax=ppR(end) + ueCoord(1) + abs(ylimnin);
        xlimmin=5*stdCoords;
        [X,Y] = meshgrid(-xlimmin:xlimmin, ylimnin:ylimnax);
        p = mvnpdf([X(:) Y(:)], ueCoord([2,1]),...
            diag([stdCoords,stdCoords].^2));
        p = reshape(p,size(X)); pcolor(X,Y,p); hold on; shading interp;
        angA = 0:360;
        for i=1:length(prbA)
            plot(ppR(i)*cosd(angA) + ueCoord(2), ...
                ppR(i)*sind(angA) + ueCoord(1), 'k-');
            text(ueCoord(2), ppR(i) + ueCoord(1), num2str(prbA(i)), ...
           'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
        plot(xPatt,yPatt,'r--','LineWidth', 2.0); hold on;
        plot(xPattDef, yPattDef, 'g-', 'LineWidth', 1.5);
        text(ueCoord(2), ueCoord(1), 'UE', 'HorizontalAlignment',...
            'center', 'VerticalAlignment', 'bottom')
        grid on; xlabel('x, м'); ylabel('y, м'); 
        axis equal; axis([-xlimmin xlimmin ylimnin ylimnax-2]);
        if ww == 0 
            legend('ДН Гаусса','Location', 'southeast');
        elseif ww == 2
            legend('Прямоугольная ДН','Location', 'southeast');
        end
        subFigNumber = subFigNumber + 1;
        ueRxPwrPlt = [ueRxPwrPlt; ueRxPwr];
    end
    sgtitle([strcat(antTypeCmt(antType), sprintf('; RMSE = %i м',s))]);
    figNumber = figNumber + 1;
end
figure; 
plot(ueRxPwrPlt(1,:).', '-'); hold on;
plot(ueRxPwrPlt(2,:).', 'x-');
plot(ueRxPwrPlt(4,:).', 'o-');
% plot(ueRxPwrPlt([1, 4, 2],:).');
grid on; xlabel('Номер точки расчета'); ylabel('P, дБ');
sgtitle([strcat(antTypeCmt(antType), sprintf('; RMSE = %i m',s))]);
legend('буз управления шириной луча',...
    'управление шириной прямоугольной ДН',...
    'управление шириной ДН Гаусса',...
    'Location', 'southeast');