clear all; close all; clc;
stdCoordsArr = [1, 3, 5, 10];  % СКО оценки координат по [x, y, z], м
NellArr = 4:1:20;              % число АЭ в одном измерении
% выбор типа антенной решетки:
% 1 - planar or uniform rectangural antenna array (URA), планарная АР
% 2 - uniform linear antenna array (ULA), линейная АР
% 3 - uniform circular antenna array (UCA); круговая АР: не поддерживается
antType = 1;
% выбор формы ДН для алгоритма управления шириной ДН (antPattCntrl = 3)
% 0 - окно Гаусса
% 1 - окно приподнятого косинуса
% 2 - прямоугольное окно
win_typeArr = [2, 0];
antTypeCmt = ["URA", "ULA", "UCA"];
for ww=1:length(win_typeArr)
    win_type = win_typeArr(ww);
    switch win_type
        case 0
            sclArr = 0.1:0.1:3.0; % окно Гаусса 
        case 1
            sclArr = 0.0:0.1:0.9; % окно приподнятого косинуса
        case 2
            sclArr = 0.1:0.1:3.0; % прямоугольное окно
    end
    c = physconst('LightSpeed');
    f = 30e9;       % несущая в диапазоне ММВ, Гц
    lamb = c/f;     % длина волны, м
    da = 0.5*c/f;   % расстояние между элементами АР
    Z = [];
    for s=stdCoordsArr
        Zn = [];
        for nel=NellArr
            Zsc = [];
            for sc=sclArr
                rng('default');
                stdCoords = s; % СКО оценки координат UE, м
                Nel = nel;     % число АЭ в одном измерении
                antElPos=createAnt(antType, Nel, da); % формирование АР
                NelFull = size(antElPos, 1);          % общее число АЭ в АР
                N = 100;  % число точек расчета
                HgNB = 0; % высота подвеса АР gNB
                Due = 50; % расстояние между gNB и UE на плоскости
                % структура параметров gNB
                gNB = createNB([0, 0, HgNB], [0, 0]);
                gNB.Steer = zeros(NelFull, 2);
                ueRxPwr = zeros(2, N);
                ueCoord = [Due, 0, 0];         % координаты UE
                gNBcoords = [gNB(:).Coords].'; % массив координат gNB
                % расстояние gNB-UE
                distSpaceT = sqrt(sum((gNB.Coords-ueCoord.').^2)); 
                BW = 2*atan2d(stdCoords,distSpaceT);% ширина ДН
                stAng = 0;                          % направление максимума
                % вектор, задающий направление из gNB в UE в глобальной СК x,y,z
                diffCoord = ueCoord - gNBcoords;
                % вектор, задающий направление из gNB в UE 
                % в локальной СК АР gNB (т.е. с учетом положения АР gNB)
                dirVect = gNB.AntOrient.'*diffCoord.';
                % расчет углов ухода от gNB к UE 
                azAng = rad2deg(atan2(dirVect(2), dirVect(1)));
                elAng = rad2deg(atan2(dirVect(3), ...
                    sqrt(sum(dirVect(1:2).^2))));
                % расчет вектора направляющих коэффициентов АР gNB
                gNB.Steer(:,1) = getAntPatternSteer(antElPos, f,...
                    azAng, elAng)/NelFull;
                % расчет вектора весовых коэфф. АР
                [w, azAngP, antPattP] = beamshapingWeight(win_type, BW,...
                                            stAng, Nel, sc);
                % применять коэффициенты w если требуемая ширина ДН не меньше, 
                % чем минимальная теоретическая ширина 0.891*lamb/Nel/da
                if (BW*sc <= rad2deg(0.891*lamb/Nel/da)...
                        || any(isnan(w)) || sum(w) == 0)
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
                for i=1:N % цикл по числу точек расчета
                    % внесение ошибки в оценку координат UE по stdCoords
                    ueCoordErr = ueCoord;
                    ueCoordErr(1:2) = ueCoordErr(1:2) +...
                        stdCoords*randn(size(ueCoord(1:2)));
                    % вектор, задающий направление из gNB в UE в глобальной 
                    % системе координат x,y,z с учетом ошибки координат UE
                    diffCoordT = ueCoordErr - gNBcoords;
                    % вектор, задающий направление из gNB в UE в системе
                    % координат АР gNB (т.е. с учетом положения АР gNB)
                    dirVectT = gNB.AntOrient.'*diffCoordT.';
                    % расчет углов отправки от gNB к UE
                    azAngT = rad2deg(atan2(dirVectT(2), dirVectT(1)));
                    elAngT = rad2deg(atan2(dirVectT(3), ...
                        sqrt(sum(dirVectT(1:2).^2))));
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
                Zi = mean(ueRxPwr(2,:) - ueRxPwr(1,:));
                fprintf('rmse = %5.1f, nel=%3i, scl=%4.2f\n', s, nel, sc);
                Zsc = [Zsc, Zi];
            end
            Zn = [Zn;Zsc];
        end
        Z{end+1} = Zn;
    end
    figure;
    [X,Y] = meshgrid(NellArr, sclArr);
    for kk=1:size(Z,2)
        subplot(2,2,kk)
        if stdCoordsArr(kk)==1
            [~,hh]=contourf(X, Y, smoothdata(Z{kk}.')); 
        else
            [~,hh]=contourf(X, Y, smoothdata(Z{kk}.'),'ShowText','on'); 
        end
        c1 = colorbar; c1.Label.String = '\DeltaP, дБ'; hold on;
        hh.LabelSpacing = 250;
        grid on; xlabel('N'); ylabel('s'); view([0, 90]); axis tight;
        title(sprintf('RMSE = %d м',stdCoordsArr(kk)));
    end
    if win_type == 0 
        sgtitle(sprintf('ДН Гаусса; %s', antTypeCmt(antType)));
    elseif win_type == 2
        sgtitle(sprintf('Прямоугольная ДН; %s', antTypeCmt(antType)));
    end   
end