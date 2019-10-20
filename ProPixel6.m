script % Projection_of_pixels_4
% Запуск всього F5
% Запуск виділеного F9

%%%%%%%%
% Дано %
%%%%%%%%
clear; clc;

% Кількість активних елементів розкладу
% в напрямку польоту (по вертикалі, Ox) TDI
pd = 5;
% поперек польоту (по горизонталі, Oy)
qd = 7;

% Номер пікселя

% Розмір чутливого елемента
% довжина [мкм]
vd = 10;
% ширина [мкм]
wd = 10;

% Період чутливих елементів
% вздовж [мкм]
Vd = 10;
% впоперек [мкм]
Wd = 10;

% Фокусна відстань об'єктива [мм]
f = 200;

% Зміщення вздовж осей координат в лінійній мірі (або в кутах)
% вздовж [град]
dAp = 0 * (pi / 180.0);
% впоперек [град]
dAq = 0 * (pi / 180.0);

% Висота орбіти  [км]
h = 250;

% Кути візування
% тангажа [град]
theta = 35 * (pi / 180.0);
% крена [град]
phi = 35 * (pi / 180.0);
% рискання [град]
psi = 0 * (pi / 180.0);

%%%%%%%%%%%%%%%%%%%%%%%%
% Уточнення параметрів %
%%%%%%%%%%%%%%%%%%%%%%%%

% Кількість опорних пікселів для відображення
% вздовж напрямку польоту
Nx0 = 3;

% поперек напрямку польоту
Ny0 = 3;

% Проміжок між пікселями в % (відсотках)
proc = 10;

% Відображення точок
% якщо id = 0 - відображати всі
% якщо id = 1 - відображати тільки вибрані
idx = 1;
idy = 1;

% Кількість точок які необхідні для опису переміщення траекторії
Nt = 10;

% Відображення графіків 1 - кожен окремо, 0 - все на одному 
graphplot = 1;

%%%%%%%%%%%%%%
% РОЗРАХУНОК %
%%%%%%%%%%%%%%

% Корекція кількості точок які необхідно відобразити
Nx = numbertoch(Nx0, pd, idx);
Ny = numbertoch(Ny0, qd, idy);

% Зміщення вздовж осей координат в лінійній мірі (або в кутах)
% вздовж [мкм]
dLp = f * tan(dAp);
% впоперек [мкм]
dLq = f * tan(dAq);

clear dAp dAq Nx0 Ny0

% Початкове лінійне значення
Lp0 = 0.5 * Vd * (1 - pd) - dLp;
Lq0 = 0.5 * Wd * (1 - qd) - dLq;

clear dLp dLq

% Визначаємо нумерацію відображуваних пікселів
Mx = masivtoch(Nx, pd, idx);
My = masivtoch(Ny, qd, idy);

% Корекція кількості точок які можливо відобразити
Nx = length(Mx);
Ny = length(My);

clear pd qd idx idy

% Розрахунок лінійних координит центрів пікселів
for i = 1:Nx
    Lpi(i) = Lp0 + (Mx(i) - 1) * Vd;
end

for i = 1:Ny
    Lqi(i) = Lq0 + (My(i) - 1) * Wd;
end

clear i Lp0 Lq0 Mx My

% Розрахунок кутових координит центрів пікселів
% Всі значення для точки "0" (0, 0)
Wxi{1} = atan(Lpi ./ f);
Wyi{1} = atan(Lqi ./ f);

clear Lpi Lqi

% Кутові величини активної частини (АЧ) пікселів
wx = anglesize(Wxi{1}, vd, f);
wy = anglesize(Wyi{1}, wd, f);

% Різниці кутових величин АЧ пікселів відносно їх центрів
dwx = diffangle(Wxi{1}, wx, vd, f);
dwy = diffangle(Wyi{1}, wy, wd, f);

clear vd wd Vd Wd f

% Розрахунок активної частини АЧ (Wxi, Wyi)
% Всі значення для точки 1 (1, -1) та (1, 0) 5
mx = 1;
Wxi{2} = koordtoch(Wxi{1}, wx, mx, dwx);
my = -1;
Wyi{2} = koordtoch(Wyi{1}, wy, my, dwy);

% Всі значення для точки 2 (1, 1) та (0, 1) 6
mx = 1;
Wxi{3} = koordtoch(Wxi{1}, wx, mx, dwx);
my = 1;
Wyi{3} = koordtoch(Wyi{1}, wy, my, dwy);

% Всі значення для точки 3 (-1, 1) та (-1, 0) 7
mx = -1;
Wxi{4} = koordtoch(Wxi{1}, wx, mx, dwx);
my = 1;
Wyi{4} = koordtoch(Wyi{1}, wy, my, dwy);

% Всі значення для точки 1 (-1, -1) та (0, -1) 8
mx = -1;
Wxi{5} = koordtoch(Wxi{1}, wx, mx, dwx);
my = -1;
Wyi{5} = koordtoch(Wyi{1}, wy, my, dwy);

clear mx my wx wy dwx dwy

% В результаті отримали масив даних для точок пікселів
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Врахуємо вплив кута рискання %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:5
    for i = 1:Nx
       for j = 1:Ny
           % Для АЧ пікселів
           theta0{k}(i, j) = atan(tan(Wxi{k}(1, i)) * cos(psi)...
               - tan(Wyi{k}(1, j)) * sin(psi));
           phi0{k}(i, j) = atan(tan(Wyi{k}(1, j)) * cos(psi)...
               + tan(Wxi{k}(1, i)) * sin(psi));
       end
    end
end

clear psi k i j Wxi Wyi

% Корегуємо кути візування згідно із врахуванням кривизни Землі
%  theta_0 = theta;
%  phi_0 = phi;
 
% Розраховуєм зміну кутових координат в залежності від кутів нахилу
for k = 1:5
    %%%%%%%%%%%%%%
    % Без нахилу %
    %%%%%%%%%%%%%%
    theta_0 = 0.0;
    phi_0 = 0.0;
    
    % спочатку по тангажу, а тоді по крену // глобальний поворот
    [theta1_0{k}, phi1_0{k}] = naklonGlob(theta_0, phi_0, theta0{k}, phi0{k}, 0);
    % спочатку по крену, а тоді по тангажу  // глобальний поворот
    [theta2_0{k}, phi2_0{k}] = naklonGlob(theta_0, phi_0, theta0{k}, phi0{k}, 1);
    % спочатку по тангажу, а тоді по крену // локальний поворот
    [theta3_0{k}, phi3_0{k}] = naklonLok(theta_0, phi_0, theta0{k}, phi0{k}, 0);
    % спочатку по крену, а тоді по тангажу // локальний поворот
    [theta4_0{k}, phi4_0{k}] = naklonLok(theta_0, phi_0, theta0{k}, phi0{k}, 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Нахил тільки по тангажу %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta_0 = atan(tan(theta) * cos(phi));
    phi_0 = 0.0;
    
    % спочатку по тангажу, а тоді по крену // глобальний поворот
    [theta1_1{k}, phi1_1{k}] = naklonGlob(theta_0, phi_0, theta0{k}, phi0{k}, 0);
    
    theta_0 = theta;
    phi_0 = 0.0;
    
    % спочатку по тангажу, а тоді по крену // локальний поворот
    [theta3_1{k}, phi3_1{k}] = naklonLok(theta_0, phi_0, theta0{k}, phi0{k}, 0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Нахил тільки по крену %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    theta_0 = 0.0;
    phi_0 = atan(tan(phi) * cos(theta));
        
    % спочатку по крену, а тоді по тангажу  // глобальний поворот
    [theta2_1{k}, phi2_1{k}] = naklonGlob(theta_0, phi_0, theta0{k}, phi0{k}, 1);
    
    theta_0 = 0.0;
    phi_0 = phi;
    
    % спочатку по крену, а тоді по тангажу // локальний поворот
    [theta4_1{k}, phi4_1{k}] = naklonLok(theta_0, phi_0, theta0{k}, phi0{k}, 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Нахил за двома кутами %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    theta_0 = theta;
    phi_0 = phi;
    
    % спочатку по тангажу, а тоді по крену // глобальний поворот
    [theta1_2{k}, phi1_2{k}] = naklonGlob(theta_0, phi_0, theta0{k}, phi0{k}, 0);
    % спочатку по крену, а тоді по тангажу  // глобальний поворот
    [theta2_2{k}, phi2_2{k}] = naklonGlob(theta_0, phi_0, theta0{k}, phi0{k}, 1);
    % спочатку по тангажу, а тоді по крену // локальний поворот
    [theta3_2{k}, phi3_2{k}] = naklonLok(theta_0, phi_0, theta0{k}, phi0{k}, 0);
    % спочатку по крену, а тоді по тангажу // локальний поворот
    [theta4_2{k}, phi4_2{k}] = naklonLok(theta_0, phi_0, theta0{k}, phi0{k}, 1);
end

% Розраховуємо траекторію переміщення центрів пікселів
for i = 1:Nx
    for j = 1:Ny
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Нахил тільки по тангажу %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        theta_0 = 0.0;
        theta_1 = atan(tan(theta) * cos(phi));
        theta_m = delta_angle(theta_0, theta_1, Nt);    % масив значень
        phi_0 = 0.0;
        
        % спочатку по тангажу, а тоді по крену // глобальний поворот
        for k = 1:Nt
            [Atheta_1{1}(i, j, k), Aphi_1{1}(i, j, k)] = naklonGlob(theta_m(k),...
                phi_0, theta0{1}(i, j), phi0{1}(i, j), 0);
        end
        
        theta_1 = theta;
        theta_m = delta_angle(theta_0, theta_1, Nt);	% масив значень
        
        % спочатку по тангажу, а тоді по крену // локальний поворот
        for k = 1:Nt
            [Atheta_3{1}(i, j, k), Aphi_3{1}(i, j, k)] = naklonLok(theta_m(k),...
                phi_0, theta0{1}(i, j), phi0{1}(i, j), 0);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Нахил тільки по крену %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        theta_0 = 0.0;
        phi_0 = 0.0;
        phi_1 = atan(tan(phi) * cos(theta));
        phi_m = delta_angle(phi_0, phi_1, Nt);          % масив значень
        
        % спочатку по крену, а тоді по тангажу  // глобальний поворот
        for k = 1:Nt
            [Atheta_2{1}(i, j, k), Aphi_2{1}(i, j, k)] = naklonGlob(theta_0,...
                phi_m(k), theta0{1}(i, j), phi0{1}(i, j), 1);
        end
        
        phi_1 = phi;
        phi_m = delta_angle(phi_0, phi_1, Nt);          % масив значень
        
        % спочатку по крену, а тоді по тангажу // локальний поворот
        for k = 1:Nt
            [Atheta_4{1}(i, j, k), Aphi_4{1}(i, j, k)] = naklonLok(theta_0,...
                phi_m(k), theta0{1}(i, j), phi0{1}(i, j), 1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Нахил за двома кутами %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        theta_0 = atan(tan(theta) * cos(phi));
        theta_1 = theta;
        theta_m = delta_angle(theta_0, theta_1, Nt);	% масив значень
        phi_0 = 0.0;
        phi_1 = phi;
        phi_m = delta_angle(phi_0, phi_1, Nt);          % масив значень
        
        % спочатку по тангажу, а тоді по крену // глобальний поворот
        for k = 1:Nt
            [Atheta_1{2}(i, j, k), Aphi_1{2}(i, j, k)] = naklonGlob(theta_m(k),...
                phi_m(k), theta0{1}(i, j), phi0{1}(i, j), 0);
        end
        
        theta_0 = 0.0;
        theta_1 = theta;
        theta_m = delta_angle(theta_0, theta_1, Nt);	% масив значень
        phi_0 = atan(tan(phi) * cos(theta));
        phi_1 = phi;
        phi_m = delta_angle(phi_0, phi_1, Nt);          % масив значень
        
        % спочатку по крену, а тоді по тангажу  // глобальний поворот
        for k = 1:Nt
            [Atheta_2{2}(i, j, k), Aphi_2{2}(i, j, k)] = naklonGlob(theta_m(k),...
                phi_m(k), theta0{1}(i, j), phi0{1}(i, j), 1);
        end
        
        theta_0 = theta;
        phi_0 = 0.0;
        phi_1 = phi;
        phi_m = delta_angle(phi_0, phi_1, Nt);          % масив значень
        
        % спочатку по тангажу, а тоді по крену // локальний поворот
        for k = 1:Nt
            [Atheta_3{2}(i, j, k), Aphi_3{2}(i, j, k)] = naklonLok(theta_0,...
                phi_m(k), theta0{1}(i, j), phi0{1}(i, j), 0);
        end
        
        theta_0 = 0.0;
        theta_1 = theta;
        theta_m = delta_angle(theta_0, theta_1, Nt);	% масив значень
        phi_0 = phi;
        
        % спочатку по крену, а тоді по тангажу // локальний поворот
        for k = 1:Nt
            [Atheta_4{2}(i, j, k), Aphi_4{2}(i, j, k)] = naklonLok(theta_m(k),...
                phi_0, theta0{1}(i, j), phi0{1}(i, j), 1);
        end
    end
end

clear i j k theta0 phi0 theta phi theta_0 phi_0 phi_1 theta_1 phi_m theta_m

% Розраховуємо лінійні координати точок
for k = 1:5
    %%%%%%%%%%%%%%
    % Без нахилу %
    %%%%%%%%%%%%%%
    Mtheta1_0{k} = h .* tan(theta1_0{k});
    Mtheta2_0{k} = h .* tan(theta2_0{k});
    Mtheta3_0{k} = h .* tan(theta3_0{k});
    Mtheta4_0{k} = h .* tan(theta4_0{k});
    
    Mphi1_0{k} = h .* tan(phi1_0{k});
    Mphi2_0{k} = h .* tan(phi2_0{k});
    Mphi3_0{k} = h .* tan(phi3_0{k});
    Mphi4_0{k} = h .* tan(phi4_0{k});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Нахил тільки по крену %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Mtheta1_1{k} = h .* tan(theta1_1{k});
    Mtheta3_1{k} = h .* tan(theta3_1{k});
    
    Mphi1_1{k} = h .* tan(phi1_1{k});
    Mphi3_1{k} = h .* tan(phi3_1{k});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Нахил тільки по крену %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Mtheta2_1{k} = h .* tan(theta2_1{k});
    Mtheta4_1{k} = h .* tan(theta4_1{k});
    
    Mphi2_1{k} = h .* tan(phi2_1{k});
    Mphi4_1{k} = h .* tan(phi4_1{k});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Нахил за двома кутами %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Mtheta1_2{k} = h .* tan(theta1_2{k});
    Mtheta2_2{k} = h .* tan(theta2_2{k});
    Mtheta3_2{k} = h .* tan(theta3_2{k});
    Mtheta4_2{k} = h .* tan(theta4_2{k});
    
    Mphi1_2{k} = h .* tan(phi1_2{k});
    Mphi2_2{k} = h .* tan(phi2_2{k});
    Mphi3_2{k} = h .* tan(phi3_2{k});
    Mphi4_2{k} = h .* tan(phi4_2{k});
end

clear k theta1_0 theta2_0 theta3_0 theta4_0 phi1_0 phi2_0 phi3_0 phi4_0...
    theta1_1 theta2_1 theta3_1 theta4_1 phi1_1 phi2_1 phi3_1 phi4_1...
    theta1_2 theta2_2 theta3_2 theta4_2 phi1_2 phi2_2 phi3_2 phi4_2...

% Розраховуємо лінійні координати траекторії
for k = 1:2
    % глобальний поворот
    Ttheta_1{k} = h .* tan(Atheta_1{k});
    Tphi_1{k} = h .* tan(Aphi_1{k});
    Ttheta_2{k} = h .* tan(Atheta_2{k});
    Tphi_2{k} = h .* tan(Aphi_2{k});
    
    % локальний поворот
    Ttheta_3{k} = h .* tan(Atheta_3{k});
    Tphi_3{k} = h .* tan(Aphi_3{k});
    Ttheta_4{k} = h .* tan(Atheta_4{k});
    Tphi_4{k} = h .* tan(Aphi_4{k});
end

clear k Atheta_1 Aphi_1 Atheta_2 Aphi_2 Atheta_3 Aphi_3 Atheta_4 Aphi_4

% Групуємо в одну змінну, щоб задати цикл перевірки на межі
for i = 1:5
    Mtheta_glob{i} = [Mtheta1_0{i}(:,:), Mtheta2_0{i}(:,:);...
        Mtheta1_1{i}(:,:), Mtheta2_1{i}(:,:);...
        Mtheta1_2{i}(:,:), Mtheta2_2{i}(:,:)];

    Mtheta_lock{i} = [Mtheta3_0{i}(:,:), Mtheta4_0{i}(:,:);...
        Mtheta3_1{i}(:,:), Mtheta4_1{i}(:,:);...
        Mtheta3_2{i}(:,:), Mtheta4_2{i}(:,:)];

    Mphi_glob{i} = [Mphi1_0{i}(:,:), Mphi2_0{i}(:,:);...
        Mphi1_1{i}(:,:), Mphi2_1{i}(:,:);...
        Mphi1_2{i}(:,:), Mphi2_2{i}(:,:)];

    Mphi_lock{i} = [Mphi3_0{i}(:,:), Mphi4_0{i}(:,:);...
        Mphi3_1{i}(:,:), Mphi4_1{i}(:,:);...
        Mphi3_2{i}(:,:), Mphi4_2{i}(:,:)];
end

for i = 2:5
    Mtheta_glob{1} = [Mtheta_glob{1}, Mtheta_glob{i}];
    Mtheta_lock{1} = [Mtheta_lock{1}, Mtheta_lock{i}];
    Mphi_glob{1} = [Mphi_glob{1}, Mphi_glob{i}];
    Mphi_lock{1} = [Mphi_lock{1}, Mphi_lock{i}];
end

Mtheta_glob = Mtheta_glob{1};
Mtheta_lock = Mtheta_lock{1};
Mphi_glob = Mphi_glob{1};
Mphi_lock = Mphi_lock{1};

if graphplot
    % Закриваємо вікна, якщо раніше було вибрано показ графіків все на одному
    close(figure(1));
    
    % Малюємо зображення проекцій
    figure(1);
    clf;
    hold on;
    grid off;   % включити сітку "on", виключити "off"
    xlabel('Поперек польоту [м]');
    ylabel('Вздовж польоту [м]');
    title('Проекція матриці пікселів (тангаж - крен) глоб. поворот');

    % Визначаємо максимум і мінімум для однакового масштабу по осям
    % глобальний поворот
    xmin = min([Mtheta_glob(:); Mtheta_lock(:)]);
    xmax = max([Mtheta_glob(:); Mtheta_lock(:)]);
    ymin = min([Mphi_glob(:); Mphi_lock(:)]);
    ymax = max([Mphi_glob(:); Mphi_lock(:)]);

    clear Mtheta_glob Mtheta_lock Mphi_glob Mphi_lock

    % Масштабуємо координатні осі
    ylen = abs(ymax - ymin);
    xlen = abs(xmax - xmin);

    ycenter = ymin + 0.5 * ylen;
    xcenter = xmin + 0.5 * xlen;

    clear ymax ymin xmax xmin

    len = max(xlen, ylen) * (1 + 0.5 * proc / 100);

    clear xlen ylen proc

    clear xlen ylen

    ylim(xcenter - 0.5 .* [len -len]);
    xlim(ycenter - 0.5 .* [len -len]);

    %  Відображаємо проекції пікселів
    for i = 1:Nx
       for j = 1:Ny
           % Без нахилу
           plot(Mphi1_0{1}(i, j), Mtheta1_0{1}(i, j), 'r+');
           plot([Mphi1_0{2}(i, j), Mphi1_0{3}(i, j),...
               Mphi1_0{4}(i, j), Mphi1_0{5}(i, j), Mphi1_0{2}(i, j)],...
               [Mtheta1_0{2}(i, j), Mtheta1_0{3}(i, j), Mtheta1_0{4}(i, j),...
               Mtheta1_0{5}(i, j), Mtheta1_0{2}(i, j)], '-k');

           % Нахил за одним кутом
           plot(Mphi1_1{1}(i, j), Mtheta1_1{1}(i, j), 'r+');
           plot([Mphi1_1{2}(i, j), Mphi1_1{3}(i, j),...
               Mphi1_1{4}(i, j), Mphi1_1{5}(i, j), Mphi1_1{2}(i, j)],...
               [Mtheta1_1{2}(i, j), Mtheta1_1{3}(i, j), Mtheta1_1{4}(i, j),...
               Mtheta1_1{5}(i, j), Mtheta1_1{2}(i, j)], '-k');

           % Нахил за двома кутами
           plot(Mphi1_2{1}(i, j), Mtheta1_2{1}(i, j), 'r+');
           plot([Mphi1_2{2}(i, j), Mphi1_2{3}(i, j),...
               Mphi1_2{4}(i, j), Mphi1_2{5}(i, j), Mphi1_2{2}(i, j)],...
               [Mtheta1_2{2}(i, j), Mtheta1_2{3}(i, j), Mtheta1_2{4}(i, j),...
               Mtheta1_2{5}(i, j), Mtheta1_2{2}(i, j)], '-k');

           % Відображаємо проекції траекторії переміщення центрів пікселів
           for t = 1:2
               for k = 1:Nt
                   Temp_t(k) = Ttheta_1{t}(i, j, k);
                   Temp_p(k) = Tphi_1{t}(i, j, k);
               end

               plot(Temp_p, Temp_t, ':b');
           end
       end
    end

    figure(2);
    clf;
    hold on;
    grid off;   % включити сітку "on", виключити "off"
    xlabel('Поперек польоту [м]');
    ylabel('Вздовж польоту [м]');
    title('Проекція матриці пікселів (крен - тангаж) глоб. поворот');

    ylim(xcenter - 0.5 .* [len -len]);
    xlim(ycenter - 0.5 .* [len -len]);

    %  Відображаємо проекції пікселів
    for i = 1:Nx
       for j = 1:Ny
           % Без нахилу
           plot(Mphi2_0{1}(i, j), Mtheta2_0{1}(i, j), 'r+');
           plot([Mphi2_0{2}(i, j), Mphi2_0{3}(i, j),...
               Mphi2_0{4}(i, j), Mphi2_0{5}(i, j), Mphi2_0{2}(i, j)],...
               [Mtheta2_0{2}(i, j), Mtheta2_0{3}(i, j), Mtheta2_0{4}(i, j),...
               Mtheta2_0{5}(i, j), Mtheta2_0{2}(i, j)], '-k');

           % Нахил за одним кутом
           plot(Mphi2_1{1}(i, j), Mtheta2_1{1}(i, j), 'r+');
           plot([Mphi2_1{2}(i, j), Mphi2_1{3}(i, j),...
               Mphi2_1{4}(i, j), Mphi2_1{5}(i, j), Mphi2_1{2}(i, j)],...
               [Mtheta2_1{2}(i, j), Mtheta2_1{3}(i, j), Mtheta2_1{4}(i, j),...
               Mtheta2_1{5}(i, j), Mtheta2_1{2}(i, j)], '-k');

           % Нахил за двома кутами
           plot(Mphi2_2{1}(i, j), Mtheta2_2{1}(i, j), 'r+');
           plot([Mphi2_2{2}(i, j), Mphi2_2{3}(i, j),...
               Mphi2_2{4}(i, j), Mphi2_2{5}(i, j), Mphi2_2{2}(i, j)],...
               [Mtheta2_2{2}(i, j), Mtheta2_2{3}(i, j), Mtheta2_2{4}(i, j),...
               Mtheta2_2{5}(i, j), Mtheta2_2{2}(i, j)], '-k');

           % Відображаємо проекції траекторії переміщення центрів пікселів
           for t = 1:2
               for k = 1:Nt
                   Temp_t(k) = Ttheta_2{t}(i, j, k);
                   Temp_p(k) = Tphi_2{t}(i, j, k);
               end

               plot(Temp_p, Temp_t, ':b');
           end
       end
    end

    figure(3);
    clf;
    hold on;
    grid off;   % включити сітку "on", виключити "off"
    xlabel('Поперек польоту [м]');
    ylabel('Вздовж польоту [м]');
    title('Проекція матриці пікселів (тангаж - крен) лок. поворот');

    ylim(xcenter - 0.5 .* [len -len]);
    xlim(ycenter - 0.5 .* [len -len]);

    %  Відображаємо проекції пікселів
    for i = 1:Nx
       for j = 1:Ny
           % Без нахилу
           plot(Mphi3_0{1}(i, j), Mtheta3_0{1}(i, j), 'r+');
           plot([Mphi3_0{2}(i, j), Mphi3_0{3}(i, j),...
               Mphi3_0{4}(i, j), Mphi3_0{5}(i, j), Mphi3_0{2}(i, j)],...
               [Mtheta3_0{2}(i, j), Mtheta3_0{3}(i, j), Mtheta3_0{4}(i, j),...
               Mtheta3_0{5}(i, j), Mtheta3_0{2}(i, j)], '-k');

           % Нахил за одним кутом
           plot(Mphi3_1{1}(i, j), Mtheta3_1{1}(i, j), 'r+');
           plot([Mphi3_1{2}(i, j), Mphi3_1{3}(i, j),...
               Mphi3_1{4}(i, j), Mphi3_1{5}(i, j), Mphi3_1{2}(i, j)],...
               [Mtheta3_1{2}(i, j), Mtheta3_1{3}(i, j), Mtheta3_1{4}(i, j),...
               Mtheta3_1{5}(i, j), Mtheta3_1{2}(i, j)], '-k');

           % Нахил за двома кутами
           plot(Mphi3_2{1}(i, j), Mtheta3_2{1}(i, j), 'r+');
           plot([Mphi3_2{2}(i, j), Mphi3_2{3}(i, j),...
               Mphi3_2{4}(i, j), Mphi3_2{5}(i, j), Mphi3_2{2}(i, j)],...
               [Mtheta3_2{2}(i, j), Mtheta3_2{3}(i, j), Mtheta3_2{4}(i, j),...
               Mtheta3_2{5}(i, j), Mtheta3_2{2}(i, j)], '-k');

           % Відображаємо проекції траекторії переміщення центрів пікселів
           for t = 1:2
               for k = 1:Nt
                   Temp_t(k) = Ttheta_3{t}(i, j, k);
                   Temp_p(k) = Tphi_3{t}(i, j, k);
               end

               plot(Temp_p, Temp_t, ':b');
           end
       end
    end

    figure(4);
    clf;
    hold on;
    grid off;   % включити сітку "on", виключити "off"
    xlabel('Поперек польоту [м]');
    ylabel('Вздовж польоту [м]');
    title('Проекція матриці пікселів (крен - тангаж) лок. поворот');

    ylim(xcenter - 0.5 .* [len -len]);
    xlim(ycenter - 0.5 .* [len -len]);

    %  Відображаємо проекції пікселів
    for i = 1:Nx
       for j = 1:Ny
           % Без нахилу
           plot(Mphi4_0{1}(i, j), Mtheta4_0{1}(i, j), 'r+');
           plot([Mphi4_0{2}(i, j), Mphi4_0{3}(i, j),...
               Mphi4_0{4}(i, j), Mphi4_0{5}(i, j), Mphi4_0{2}(i, j)],...
               [Mtheta4_0{2}(i, j), Mtheta4_0{3}(i, j), Mtheta4_0{4}(i, j),...
               Mtheta4_0{5}(i, j), Mtheta4_0{2}(i, j)], '-k');

           % Нахил за одним кутом
           plot(Mphi4_1{1}(i, j), Mtheta4_1{1}(i, j), 'r+');
           plot([Mphi4_1{2}(i, j), Mphi4_1{3}(i, j),...
               Mphi4_1{4}(i, j), Mphi4_1{5}(i, j), Mphi4_1{2}(i, j)],...
               [Mtheta4_1{2}(i, j), Mtheta4_1{3}(i, j), Mtheta4_1{4}(i, j),...
               Mtheta4_1{5}(i, j), Mtheta4_1{2}(i, j)], '-k');

           % Нахил за двома кутами
           plot(Mphi4_2{1}(i, j), Mtheta4_2{1}(i, j), 'r+');
           plot([Mphi4_2{2}(i, j), Mphi4_2{3}(i, j),...
               Mphi4_2{4}(i, j), Mphi4_2{5}(i, j), Mphi4_2{2}(i, j)],...
               [Mtheta4_2{2}(i, j), Mtheta4_2{3}(i, j), Mtheta4_2{4}(i, j),...
               Mtheta4_2{5}(i, j), Mtheta4_2{2}(i, j)], '-k');

           % Відображаємо проекції траекторії переміщення центрів пікселів
           for t = 1:2
               for k = 1:Nt
                   Temp_t(k) = Ttheta_4{t}(i, j, k);
                   Temp_p(k) = Tphi_4{t}(i, j, k);
               end

               plot(Temp_p, Temp_t, ':b');
           end
       end
    end
    
else
    % Закриваємо вікна, якщо раніше було вибрано показ графіків окремо
    close(figure(1), figure(2), figure(3), figure(4));
    
    % Малюємо зображення проекцій
    figure(1);
    clf;
    
    % Розбиваємо графіки на 2x2 поля, поле №1
    subplot(2, 2, 1);
    
    hold on;
    grid off;   % включити сітку "on", виключити "off"
    xlabel('Поперек польоту [м]');
    ylabel('Вздовж польоту [м]');
    title('(тангаж - крен) глоб. поворот');

    % Визначаємо максимум і мінімум для однакового масштабу по осям
    % глобальний поворот
    xmin = min([Mtheta_glob(:); Mtheta_lock(:)]);
    xmax = max([Mtheta_glob(:); Mtheta_lock(:)]);
    ymin = min([Mphi_glob(:); Mphi_lock(:)]);
    ymax = max([Mphi_glob(:); Mphi_lock(:)]);

    clear Mtheta_glob Mtheta_lock Mphi_glob Mphi_lock

    % Масштабуємо координатні осі
    ylen = abs(ymax - ymin);
    xlen = abs(xmax - xmin);

    ycenter = ymin + 0.5 * ylen;
    xcenter = xmin + 0.5 * xlen;

    clear ymax ymin xmax xmin

    len = max(xlen, ylen) * (1 + 0.5 * proc / 100);

    clear xlen ylen proc

    clear xlen ylen

    ylim(xcenter - 0.5 .* [len -len]);
    xlim(ycenter - 0.5 .* [len -len]);

    %  Відображаємо проекції пікселів
    for i = 1:Nx
       for j = 1:Ny
           % Без нахилу
           plot(Mphi1_0{1}(i, j), Mtheta1_0{1}(i, j), 'r+');
           plot([Mphi1_0{2}(i, j), Mphi1_0{3}(i, j),...
               Mphi1_0{4}(i, j), Mphi1_0{5}(i, j), Mphi1_0{2}(i, j)],...
               [Mtheta1_0{2}(i, j), Mtheta1_0{3}(i, j), Mtheta1_0{4}(i, j),...
               Mtheta1_0{5}(i, j), Mtheta1_0{2}(i, j)], '-k');

           % Нахил за одним кутом
           plot(Mphi1_1{1}(i, j), Mtheta1_1{1}(i, j), 'r+');
           plot([Mphi1_1{2}(i, j), Mphi1_1{3}(i, j),...
               Mphi1_1{4}(i, j), Mphi1_1{5}(i, j), Mphi1_1{2}(i, j)],...
               [Mtheta1_1{2}(i, j), Mtheta1_1{3}(i, j), Mtheta1_1{4}(i, j),...
               Mtheta1_1{5}(i, j), Mtheta1_1{2}(i, j)], '-k');

           % Нахил за двома кутами
           plot(Mphi1_2{1}(i, j), Mtheta1_2{1}(i, j), 'r+');
           plot([Mphi1_2{2}(i, j), Mphi1_2{3}(i, j),...
               Mphi1_2{4}(i, j), Mphi1_2{5}(i, j), Mphi1_2{2}(i, j)],...
               [Mtheta1_2{2}(i, j), Mtheta1_2{3}(i, j), Mtheta1_2{4}(i, j),...
               Mtheta1_2{5}(i, j), Mtheta1_2{2}(i, j)], '-k');

           % Відображаємо проекції траекторії переміщення центрів пікселів
           for t = 1:2
               for k = 1:Nt
                   Temp_t(k) = Ttheta_1{t}(i, j, k);
                   Temp_p(k) = Tphi_1{t}(i, j, k);
               end

               plot(Temp_p, Temp_t, ':b');
           end
       end
    end

    % Розбиваємо графіки на 2x2 поля, поле №2
    subplot(2, 2, 2);
    
    hold on;
    grid off;   % включити сітку "on", виключити "off"
    xlabel('Поперек польоту [м]');
    ylabel('Вздовж польоту [м]');
    title('(крен - тангаж) глоб. поворот');

    ylim(xcenter - 0.5 .* [len -len]);
    xlim(ycenter - 0.5 .* [len -len]);

    %  Відображаємо проекції пікселів
    for i = 1:Nx
       for j = 1:Ny
           % Без нахилу
           plot(Mphi2_0{1}(i, j), Mtheta2_0{1}(i, j), 'r+');
           plot([Mphi2_0{2}(i, j), Mphi2_0{3}(i, j),...
               Mphi2_0{4}(i, j), Mphi2_0{5}(i, j), Mphi2_0{2}(i, j)],...
               [Mtheta2_0{2}(i, j), Mtheta2_0{3}(i, j), Mtheta2_0{4}(i, j),...
               Mtheta2_0{5}(i, j), Mtheta2_0{2}(i, j)], '-k');

           % Нахил за одним кутом
           plot(Mphi2_1{1}(i, j), Mtheta2_1{1}(i, j), 'r+');
           plot([Mphi2_1{2}(i, j), Mphi2_1{3}(i, j),...
               Mphi2_1{4}(i, j), Mphi2_1{5}(i, j), Mphi2_1{2}(i, j)],...
               [Mtheta2_1{2}(i, j), Mtheta2_1{3}(i, j), Mtheta2_1{4}(i, j),...
               Mtheta2_1{5}(i, j), Mtheta2_1{2}(i, j)], '-k');

           % Нахил за двома кутами
           plot(Mphi2_2{1}(i, j), Mtheta2_2{1}(i, j), 'r+');
           plot([Mphi2_2{2}(i, j), Mphi2_2{3}(i, j),...
               Mphi2_2{4}(i, j), Mphi2_2{5}(i, j), Mphi2_2{2}(i, j)],...
               [Mtheta2_2{2}(i, j), Mtheta2_2{3}(i, j), Mtheta2_2{4}(i, j),...
               Mtheta2_2{5}(i, j), Mtheta2_2{2}(i, j)], '-k');

           % Відображаємо проекції траекторії переміщення центрів пікселів
           for t = 1:2
               for k = 1:Nt
                   Temp_t(k) = Ttheta_2{t}(i, j, k);
                   Temp_p(k) = Tphi_2{t}(i, j, k);
               end

               plot(Temp_p, Temp_t, ':b');
           end
       end
    end

    % Розбиваємо графіки на 2x2 поля, поле №3
    subplot(2, 2, 3);
    
    hold on;
    grid off;   % включити сітку "on", виключити "off"
    xlabel('Поперек польоту [м]');
    ylabel('Вздовж польоту [м]');
    title('(тангаж - крен) лок. поворот');

    ylim(xcenter - 0.5 .* [len -len]);
    xlim(ycenter - 0.5 .* [len -len]);

    %  Відображаємо проекції пікселів
    for i = 1:Nx
       for j = 1:Ny
           % Без нахилу
           plot(Mphi3_0{1}(i, j), Mtheta3_0{1}(i, j), 'r+');
           plot([Mphi3_0{2}(i, j), Mphi3_0{3}(i, j),...
               Mphi3_0{4}(i, j), Mphi3_0{5}(i, j), Mphi3_0{2}(i, j)],...
               [Mtheta3_0{2}(i, j), Mtheta3_0{3}(i, j), Mtheta3_0{4}(i, j),...
               Mtheta3_0{5}(i, j), Mtheta3_0{2}(i, j)], '-k');

           % Нахил за одним кутом
           plot(Mphi3_1{1}(i, j), Mtheta3_1{1}(i, j), 'r+');
           plot([Mphi3_1{2}(i, j), Mphi3_1{3}(i, j),...
               Mphi3_1{4}(i, j), Mphi3_1{5}(i, j), Mphi3_1{2}(i, j)],...
               [Mtheta3_1{2}(i, j), Mtheta3_1{3}(i, j), Mtheta3_1{4}(i, j),...
               Mtheta3_1{5}(i, j), Mtheta3_1{2}(i, j)], '-k');

           % Нахил за двома кутами
           plot(Mphi3_2{1}(i, j), Mtheta3_2{1}(i, j), 'r+');
           plot([Mphi3_2{2}(i, j), Mphi3_2{3}(i, j),...
               Mphi3_2{4}(i, j), Mphi3_2{5}(i, j), Mphi3_2{2}(i, j)],...
               [Mtheta3_2{2}(i, j), Mtheta3_2{3}(i, j), Mtheta3_2{4}(i, j),...
               Mtheta3_2{5}(i, j), Mtheta3_2{2}(i, j)], '-k');

           % Відображаємо проекції траекторії переміщення центрів пікселів
           for t = 1:2
               for k = 1:Nt
                   Temp_t(k) = Ttheta_3{t}(i, j, k);
                   Temp_p(k) = Tphi_3{t}(i, j, k);
               end

               plot(Temp_p, Temp_t, ':b');
           end
       end
    end

    % Розбиваємо графіки на 2x2 поля, поле №4
    subplot(2, 2, 4);
    
    hold on;
    grid off;   % включити сітку "on", виключити "off"
    xlabel('Поперек польоту [м]');
    ylabel('Вздовж польоту [м]');
    title('(крен - тангаж) лок. поворот');

    ylim(xcenter - 0.5 .* [len -len]);
    xlim(ycenter - 0.5 .* [len -len]);

    %  Відображаємо проекції пікселів
    for i = 1:Nx
       for j = 1:Ny
           % Без нахилу
           plot(Mphi4_0{1}(i, j), Mtheta4_0{1}(i, j), 'r+');
           plot([Mphi4_0{2}(i, j), Mphi4_0{3}(i, j),...
               Mphi4_0{4}(i, j), Mphi4_0{5}(i, j), Mphi4_0{2}(i, j)],...
               [Mtheta4_0{2}(i, j), Mtheta4_0{3}(i, j), Mtheta4_0{4}(i, j),...
               Mtheta4_0{5}(i, j), Mtheta4_0{2}(i, j)], '-k');

           % Нахил за одним кутом
           plot(Mphi4_1{1}(i, j), Mtheta4_1{1}(i, j), 'r+');
           plot([Mphi4_1{2}(i, j), Mphi4_1{3}(i, j),...
               Mphi4_1{4}(i, j), Mphi4_1{5}(i, j), Mphi4_1{2}(i, j)],...
               [Mtheta4_1{2}(i, j), Mtheta4_1{3}(i, j), Mtheta4_1{4}(i, j),...
               Mtheta4_1{5}(i, j), Mtheta4_1{2}(i, j)], '-k');

           % Нахил за двома кутами
           plot(Mphi4_2{1}(i, j), Mtheta4_2{1}(i, j), 'r+');
           plot([Mphi4_2{2}(i, j), Mphi4_2{3}(i, j),...
               Mphi4_2{4}(i, j), Mphi4_2{5}(i, j), Mphi4_2{2}(i, j)],...
               [Mtheta4_2{2}(i, j), Mtheta4_2{3}(i, j), Mtheta4_2{4}(i, j),...
               Mtheta4_2{5}(i, j), Mtheta4_2{2}(i, j)], '-k');

           % Відображаємо проекції траекторії переміщення центрів пікселів
           for t = 1:2
               for k = 1:Nt
                   Temp_t(k) = Ttheta_4{t}(i, j, k);
                   Temp_p(k) = Tphi_4{t}(i, j, k);
               end

               plot(Temp_p, Temp_t, ':b');
           end
       end
    end
end

clc;

clear all
