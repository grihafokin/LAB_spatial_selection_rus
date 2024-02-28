% функция расчета вектора весовых коэффициентов для синтеза диаграммы 
% направленности (ДН) заданной формы и ширины главного лепестка луча
% входные параметры: 
% win_type - тип оконной функции для синтеза ДН 
%            (прямоугольной, приподнятый косинус, функция Гаусса)
% BW -       ширина главно лепестка ДН, градусы
% stAng -    направление максимума ДН в азимутальной плоскости, градусы
% N -        число АЭ
% scl -      множитель ширины главного лепестка ДН
% выходные параметры:
% w -        вектор весовых коэффициентов
% antPattP - заданная форма ДН
% azAngP - значения азимутальных углов, в которых задана antPattP
function [w, azAngP, antPattP]  = ...
    beamshapingWeight(win_type, BW, stAng, N, scl, varargin)
k = 0:N-1;
% расчет точек, в которых будет задаваться синтезируемая ДН, формула (34)
psi_k = wrapToPi((k - (N - 1)/2)*2*pi/N);
% пересчет psi_k в азимутальные углы; для d = lambda/2 
% формула (29) упрощается 2*pi/lamb*sin(az)*d = pi*sin(az); 
% знак "-" для удобства представления данных
azAngW = -asind(psi_k/pi);
Na = length(azAngW);
azAngP = -90:0.1:89;
Np = length(azAngP);
antPattP = zeros(1, Np);
% формирование ДН заданной формы и ширины
switch win_type
    case 0  % ДН в виде функции Гаусса (см. формулу (43))
        antPattW = exp(-((azAngW-stAng).^2/2/(BW*scl/2)^2));
        antPattP = exp(-((azAngP-stAng).^2/2/(BW*scl/2)^2));
    case 1 % ДН в виде функции приподнятого косинуса 
           % (в моделировании не используется)
        beta = 0.35;
        B = BW*scl/(1-beta);
        antPattW = zeros(1, Na);
        for j=1:Na
            if (azAngW(j) > -(1-beta)*B/2+stAng && ...
                    azAngW(j) < (1-beta)*B/2+stAng)
                antPattW(j) = 1;
            elseif (azAngW(j) > (1-beta)*B/2+stAng &&...
                    azAngW(j) < (1+beta)*B/2+stAng) ||...
                    (azAngW(j) < -(1-beta)*B/2+stAng &&...
                    azAngW(j) > -(1+beta)*B/2+stAng)
                antPattW(j) = 0.5*(1+cos(pi/beta/B*(abs(azAngW(j)-stAng)...
                    - (1-beta)*B/2)));
            end
        end
        for j=1:Np
            if (azAngP(j) > -(1-beta)*B/2+stAng &&...
                    azAngP(j) < (1-beta)*B/2+stAng)
                antPattP(j) = 1;
            elseif (azAngP(j) > (1-beta)*B/2+stAng &&...
                    azAngP(j) < (1+beta)*B/2+stAng) ||...
                    (azAngP(j) < -(1-beta)*B/2+stAng &&...
                    azAngP(j) > -(1+beta)*B/2+stAng)
                antPattP(j) = 0.5*(1+cos(pi/beta/B*(abs(azAngP(j)-stAng)...
                    - (1-beta)*B/2)));
            end
        end
    case 2 % прямоугольная ДН (формула (42))
        antPattW = zeros(size(azAngW));
        antPattW(azAngW>(-BW*scl/2+stAng) & azAngW<(BW*scl/2+stAng)) = 1;
        % дополнительная обработка случая, когда в antPattW только одно
        % не нулевое значение; в этом случае соседние точки слева и
        % справа интерполируются линейно
        rEdge = (-BW*scl/2+stAng);
        lEdge = (BW*scl/2+stAng);
        if (sum(antPattW) == 1)
            lInd = find(azAngW>lEdge,1,'last');
            rInd = find(azAngW<rEdge,1,'first');
            if (lInd == 1)
                antPattW(lInd) = 1 - (azAngW(lInd) - lEdge)/...
                    (azAngW(lInd) - azAngW(lInd+1))/2;
            elseif (lInd > 1)
                antPattW(lInd) = 1 - (azAngW(lInd) - lEdge)/...
                    (azAngW(lInd-1) - lEdge);
            end
            if (rInd == N)
                antPattW(rInd) = 1 - (azAngW(rInd) - rEdge)/...
                    (azAngW(rInd) - azAngW(rInd-1))/2;
            elseif (rInd < N)
                antPattW(rInd) = 1 - (azAngW(rInd) - rEdge)/...
                    (azAngW(rInd+1) - rEdge);
            end
        end
        antPattP(azAngP>(-BW*scl/2+stAng) & azAngP<(BW*scl/2+stAng)) = 1;
end 
Ak = conj(antPattW).*exp(-1i*psi_k*(N-1)/2); % расчет Ak (формула (37))
n = 0:N-1;
% расчет bn (формула (36)); для ускорения работы используется ifft
bn = ifft(Ak); % bn = Ak*exp(1i*(k'*n)*2*pi/N)/N;
% расчет вектора коэфф. (формула (40))
w = (bn.*exp(-1i*n*pi*(N-1)/N))';
% дополнительная обработка вектора коэффициентов w 
% оконной функцией Хэмминга (формула (44))
if (nargin == 6)
    useW = varargin{1};
else
    useW = 1;
end
if win_type == 2 && useW == 1
    alph = 0.54;
    % коэффициенты окна Хэмминга
    wW = alph - (1-alph)*cos(2*pi*(0:N-1)/(N-1));
    w = w.*wW.';
end
end