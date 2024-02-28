% функция создания сценария расчета (положение gNB, траектории движения UE)
function [gNB, ueNode, d, T, v] = createScenarion(sceneN, varargin)
T = 0.1;  % период измерений
v = 10;   % скорость движения устройства UE, м/с
switch sceneN
    % gNB расположены на одной стороне относительно траектории
    % движения UE. UE двигаются параллельно на расстоянии d
    case 1 
        % расстояние между траекториями движения двух UE;
        % в данном сценарии траектории UE разнесены на d по Y
        if (nargin == 2)
            d = varargin{1};
        else
            d = 0;
        end
        % создание массива из двух структур gNB; параметры gNB в createNB
        gNB(1) = createNB([25, 50, 5], [-90, -1]);
        gNB(2) = createNB([125, 50, 5], [-90, -1]);
        % создание массива из двух структур UE; параметры UE в createUEnode
        ueNode(1) = createUEnode([0; 150], [0; 0], 0, 1, v, T, [90, 0]);
        ueNode(2) = createUEnode([150; 0], [0; 0], 0, 2, v, T, [90, 0]); 
         % построения набора траекторий движения UE для разных значений d; 
         % траектория UE2 не меняется, на d смещается траектория UE1
         Nd = length(d); % число значений d
         trajUE1 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
         trajUE2 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
         for i=1:Nd
             trajUE1(:,:,i) = ueNode(1).Trajectory + [0, d(i), 0];
             trajUE2(:,:,i) = ueNode(2).Trajectory;
         end
         % компоновка набора траекторий в один общий массив координат
         % для каждой UE
         trajUE1 = permute(trajUE1, [1 3 2]);
         trajUE1 = reshape(trajUE1, [],3 ,1);
         trajUE2 = permute(trajUE2, [1 3 2]);
         trajUE2 = reshape(trajUE2, [],3 ,1);
         % сохранение траектории в структуре UE
         ueNode(1).Trajectory = trajUE1;
         ueNode(2).Trajectory = trajUE2;
            
    % gNB расположены по разные стороны относительно траектории
    % движения UE; UE двигаются параллельно на расстоянии d
    case 2
        % расстояние между траекториями движения двух UE;
        % в данном сценарии траектории UE разнесены на d по Y
        d = 0:1:10;
        % создание массива из двух структур gNB
        gNB(1) = createNB([25, 50, 5], [-90, -1]);
        gNB(2) = createNB([125, 0, 5], [90, -1]);
        % создание массива из двух структур UE
        ueNode(1) = createUEnode([0; 150], [20; 20], 0, 1, v, T, [90, 0]);
        ueNode(2) = createUEnode([150; 0], [20; 20], 0, 2, v, T, [90, 0]);
        % построения набора траекторий движения UE для разных значений d;
        % траектория UE2 не меняется, на d смещается траектория UE1
        Nd = length(d); % число значений d
        trajUE1 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        trajUE2 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        for i=1:Nd
            trajUE1(:,:,i) = ueNode(1).Trajectory + [0, d(i), 0];
            trajUE2(:,:,i) = ueNode(2).Trajectory;
        end
        % компоновка набора траекторий в один общий массив координат
        % для каждой UE
        trajUE1 = permute(trajUE1, [1 3 2]);
        trajUE1 = reshape(trajUE1, [],3 ,1);
        trajUE2 = permute(trajUE2, [1 3 2]);
        trajUE2 = reshape(trajUE2, [],3 ,1);
        % сохранение траектории в структуре UE
        ueNode(1).Trajectory = trajUE1;
        ueNode(2).Trajectory = trajUE2;
    % gNB расположены на одной стороне относительно траектории
    % движения UE; UE двигаются друг за другом на расстоянии d
    case 3
        % расстояние между траекториями движения двух UE;
        % в данном сценарии траектории UE разнесены на d по X
        d = 0:2:150;
        % создание массива из двух структур gNB
        gNB(1) = createNB([25, 50, 5], [-90, -1]);
        gNB(2) = createNB([125, 50, 5], [-90, -1]);
        % создание массива из двух структур UE
        ueNode(1) = createUEnode([0; 150], [0; 0], 0, 1, v, T, [90, 0]);
        ueNode(2) = createUEnode([0; 150], [0; 0], 0, 2, v, T, [90, 0]);
            
        % построения набора траекторий движения UE для разных значений d;
        % траектория UE2 не меняется, на d смещается траектория UE1
        Nd = length(d); % число значений d
        trajUE1 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        trajUE2 = zeros(length(ueNode(1).Trajectory(:,1)), 3, Nd);
        for i=1:Nd
            trajUE1(:,:,i) = ueNode(1).Trajectory + [d(i), (i-1)/1e10, 0];
            trajUE2(:,:,i) = ueNode(2).Trajectory;
        end
        % компоновка набора траекторий в один общий массив координат
        % для каждой UE
        trajUE1 = permute(trajUE1, [1 3 2]);
        trajUE1 = reshape(trajUE1, [],3 ,1);
        trajUE2 = permute(trajUE2, [1 3 2]);
        trajUE2 = reshape(trajUE2, [],3 ,1);
        % сохранение траектории в структуре UE
        ueNode(1).Trajectory = trajUE1;
        ueNode(2).Trajectory = trajUE2;
    otherwise
end
end