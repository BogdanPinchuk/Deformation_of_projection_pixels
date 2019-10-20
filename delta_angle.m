function masiv = delta_angle(angle0, angle1, Nt)
% Функція розрахунку нахилу КА

% angle0 - початковий кут перед відхиленням
% angle1 - кінечений кут після відхилення
% Nt - кількість точок на які необхідно розбити траекторію
% masiv - масив розрахованих точок

% Перевірка на кількість точок
if Nt < 2
    Nt = 2;
end

% Визначаємо крок
delta = (angle1 - angle0) ./ (Nt - 1);

for i = 1:Nt
    masiv(i) = angle0 + delta .* (i - 1);
end

end
