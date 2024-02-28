clear all; close all; clc;
c = physconst('LightSpeed');
f = 30e9;       % несущая в диапазоне ММВ, Гц
lamb = c/f;     % длина волны, м
da = 0.5*c/f;   % расстояние между элементами АР
BW = 50;        % ширина луча HPBW, градусы
Nel = 20;       % число АЭ в одном измерении
% selection of antenna array type
% 1 - planar or uniform rectangural antenna array (URA), планарная АР
% 2 - uniform linear antenna array (ULA), линейная АР
% 3 - uniform circular antenna array (UCA); круговая АР, не поддерживается
antType = 1;
antElPos = createAnt(antType, Nel, da); % формирование АР

figNumber = 1;
for j=[0,1]
    % расчет коэфф. для прямоугольной ДН
    [wr, azAngP, antPattPr] = beamshapingWeight(2, BW, 0, Nel, 1, j);
    if (antType == 1)
        % расчет вектора весовых коэфф. планарной АР для вертикальных АЭ
        wr = repmat(wr, Nel, 1)/Nel;
        wr = wr(:);
    end
    gr = zeros(1, length(azAngP));
    for i=1:length(azAngP)
        gr(i) = getAntPatternG(antElPos, f, azAngP(i), 0, wr, 0);
    end

    figure;
    plot(azAngP, antPattPr,'LineWidth', 2); hold on; 
    plot(azAngP, gr,'--','LineWidth', 2); grid on;
    xlabel('\phi, \circ'); ylabel('|A(\phi)|');
    legend('задано', 'синтезировано');
    figNumber = figNumber + 1;
end
% расчет коэфф. для ДН Гаусса
[wg, azAngP, antPattPg] = beamshapingWeight(0, BW, 0, Nel, 1);
if (antType == 1)
    % расчет вектора весовых коэфф. планарной АР для вертикальных АЭ
    wg = repmat(wg, Nel, 1)/Nel;
    wg = wg(:);
end
gg = zeros(1, length(azAngP));
for i=1:length(azAngP)
    gg(i) = getAntPatternG(antElPos, f, azAngP(i), 0, wg, 0);
end

figure;
plot(azAngP, antPattPg, 'LineWidth', 2); hold on;
plot(azAngP, gg, '--', 'LineWidth', 2); grid on;
xlabel('\phi, \circ'); ylabel('|A(\phi)|');
legend('задано', 'синтезировано');