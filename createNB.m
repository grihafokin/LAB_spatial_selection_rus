% функция создания структуры gNB
% входные параметры:
% coords - координаты [x,y,z] gNB, м
% antDir - ориентация АР gNB (направление оси симметрии) 
%          [азимут, наклон], градусы
function gNB = createNB(coords, antDir)
gNB.Coords = coords.'; % координаты gNB
gNB.AntDir = antDir.'; % ориентация АР gNB
% матрица поворота согласно значениям antDir, используется для пересчета 
% векторов направлений из глобальных координат в систему координат АР gNB
gNB.AntOrient = rotz(antDir(1))*roty(-antDir(2));
gNB.Steer = 1;
gNB.UEPolarCoord = [0, 0, 0];
end