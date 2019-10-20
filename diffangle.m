function dw = diffangle(W, w, s, f)
% ������� ���������� ������ ������� ��������� ������

% dw - ������ ������� ��������� ������
% w - ������ �������� ������� ������� ������ ��� ������ ������
% W - Wxi, Wyi - ����� ����� ������ ������
% s - vd, wd, Vd, Wd - ����� ������� ������� ������ ��� ������ ������
% f - ������� �������

dw = sign(W) .* acos(((s ./ (2 .* f)) .* (cos(W) .^ 2)...
    .* sin(w) + cos(w)));

end
