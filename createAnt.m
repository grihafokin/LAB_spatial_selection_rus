% функция создания антенной решетки 
% входные параметры:
% antType - тип АР
% Nel     - число АЭ
% d       - расстояние между АЭ, м
% выходные параметры:
% antPos  - массив координат [x,y,z] АЭ, м 
% (размерность [Nel x 3] для линейной АР, [Nel^2 x 3] для планарной АР
function antPos = createAnt(antType, Nel, d)
switch antType
    case 1 % планарная (квадратная) АР
        % число антенных элементов на одной стороне квадратной 
        % антенной решетки; общее число элементов Nel^2
        yCoords = repmat((-(Nel-1)/2:1:(Nel-1)/2).'*d, Nel, 1);
        zCoords = reshape(repmat((-(Nel-1)/2:1:(Nel-1)/2)*d, Nel, 1),[],1);
        zCoords = zCoords(:);
        xCoords = zeros(size(yCoords));
        antPos = [xCoords, yCoords, zCoords];        
    case 2 % линейная АР
        % линейная антенная решетка из Nel элементов
        antLen = (Nel-1)*d;
        yCoords = (-antLen/2:d:antLen/2).';
        zCoords = zeros(size(yCoords));
        xCoords = zeros(size(yCoords));
        antPos = [xCoords, yCoords, zCoords];  
    case 3 % круговая АР  
        % общее число элементов Nel = Nel*2;
        % расчет радиуса круговой решетки через формулу длины хорды дуги 
        % окружности; длина хорды равна заданному расстоянию между АЭ d
        r = d/2/sind(360/Nel/2);
        antDph = (0:360/Nel:360 - 360/Nel).';
        xCoords = r*cosd(antDph);
        yCoords = r*sind(antDph);
        zCoords = zeros(size(yCoords));
        antPos = [xCoords, yCoords, zCoords];  
end
end