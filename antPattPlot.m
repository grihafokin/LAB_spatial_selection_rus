% функция построения ДН АР
% antElPos   - массив координат [x,y,z] АЭ, м
% f          - рабочая частота, Гц
% nodeStruct - структура с параметрами узла (gNB или UE)
% angStep    - шаг сетки углов по азимуту и углу места
% scl        - коэфф. масштабирования ДН (для удобства визуализации)
% backLobe   - использовать подавление обратного лепестка
% plt2D      - отображать ДН в 2D на срезе угла места равном 0
% [x, y, z]  - координаты поверхности для отображения ДН (безразмерная 
%                        величина, используется только для отображения)
function [x, y, z, pp] = antPattPlot(antElPos, f, nodeStruct, ...
                                     angStep, scl, backLobe, plt2D)   
if (nargin == 6)
    plt2D = 0;
end
    
% сетка углов по азимуту и углу места для расчета ДН
azA = 0:angStep:360;
if (plt2D == 0)
    elA = -90:angStep:90;
else
    elA = 0;
end
% инициализация массивов для сохранения значений ДН
azN = length(azA);
elN = length(elA);
x = zeros(elN, azN);
y = zeros(elN, azN);
z = zeros(elN, azN);
pp = zeros(elN, azN);
    
for i=1:elN % цикл по массиву углов места
    for j=1:azN % цикл по массиву азимутов
        % расчет значения ДН АР для i-го угла места и j-го азимута 
        % с учетом вектора направляющих коэфф. АР (nodeStruct.Steer)
        if (plt2D == 0)
            p = getAntPatternG(antElPos, f, azA(j), elA(i),...
                nodeStruct.Steer, backLobe)*scl;
        else
            p = getAntPatternG(antElPos, f, azA(j),...
                nodeStruct.UEPolarCoord(2), ...
                nodeStruct.Steer, backLobe)*scl;
        end
        pp(i,j) = p;
        
        % пересчет значений ДН АР из полярных координат 
        % в прямоугольные с учетом ориентации АР (nodeStruct.AntOrient)
        r = p.*cosd(elA(i));
        if (plt2D == 0)
            xyz = nodeStruct.AntOrient*[r.*cosd(azA(j)); ...
                r.*sind(azA(j)); p.*sind(elA(i))];
        else
          xyz = [[cosd(-nodeStruct.AntDir(1)),sind(-nodeStruct.AntDir(1)); 
             -sind(-nodeStruct.AntDir(1)), cosd(-nodeStruct.AntDir(1))]*...
             [r.*cosd(azA(j)); r.*sind(azA(j))];0];
        end
        x(i,j) = xyz(1);
        y(i,j) = xyz(2);
        z(i,j) = xyz(3);
    end % for j=1:azN % цикл по массиву азимутов
end % for i=1:elN % цикл по массиву углов места
end