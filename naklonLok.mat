function [Theta, Phi] = naklonLok(theta, phi, Wx, Wy, ind)
% Функція розрахунку нахилу КА

% Theta - theta1
% Phi - phi1
% theta, phi - стандарт
% Wx - theta0
% Wy - phi0
% ind - індекс який визначає як рухатиметься система
% якщо 0 - то спочатку по тангажу, а тоді по крену
% якщо 1 - то спочатку по крену, а тоді по тангажу

% Попередній розрахунок
theta_2 = atan(tan(phi) .* cos(theta));
phi_2 = atan(tan(theta) .* cos(phi));

Wx_2 = atan(tan(Wy) .* cos(Wx));
Wy_2 = atan(tan(Wx) .* cos(Wy));

if ind == 0
   % якщо ind = 0
   % для локального повороту
   dtheta = atan(tan(Wy_2) ./ cos(theta_2 + Wy));
   Theta = theta + dtheta;
   Phi = atan(tan(theta_2 + Wy) .* cos(dtheta) ./ cos(Theta));
else
   % якщо ind = 1
   % для локального повороту
   dphi = atan(tan(Wx_2) ./ cos(phi_2 + Wx));
   Phi = phi + dphi;
   Theta = atan(tan(phi_2 + Wx) .* cos(dphi) ./ cos(Phi));
end

Theta = real(Theta);
Phi = real(Phi);

end
