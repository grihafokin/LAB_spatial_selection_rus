% функция построения траектории движения UE
% входные параметры:
% trajPnts - начальныя и конечные значения траектории UE по X, Y, Z, м
% v - скорость движения UE, м/с
% T - период измерения координат UE, с
function xyUN = getTrajectory(trajPnts, v, T)
% общее время сценария движения UE
times = [0; sqrt(sum((trajPnts(2:end, :) - trajPnts(1:end-1, :)).^2,2))/v];
% время в каждой опорной точке (trajPnts), описывающих траекторию движения
elapsedTime = zeros(1, length(times));    
for i=1:length(times)
    elapsedTime(i) = sum(times(1:i));
end
% встроенные функции для формирования координат 
% траектории по заданным опорным точкам
ts = trackingScenario('UpdateRate', 1/T);
target = platform(ts);
traj = waypointTrajectory('Waypoints', trajPnts, ...
    'TimeOfArrival', elapsedTime, 'SampleRate', 1/T);
target.Trajectory = traj;
r = record(ts);
posUN = [r(:).Poses];
% координаты точек траектории UE
xyUN = vertcat(posUN.Position);
end