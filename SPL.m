%L0325042_RafiNazwanMaulana
clc;
clear;
format long;

% Sistem persamaan linear
A = [ 5 -1 -1  0  0;
     -1  5 -1 -1  0;
     -1 -1  4 -1 -1;
      0  0  1  4 -2;
      0  1 -1  1  4];

b = [-1; 2; 6; 2; -1];

n = length(b);
maks_iterasi = 20;
toleransi = 0.001;

% METODE GAUSS
Ag = A;
bg = b;
xGauss = zeros(n,1);

tic;
for k = 1:n-1
    for i = k+1:n
        faktor = Ag(i,k) / Ag(k,k);
        Ag(i,k:n) = Ag(i,k:n) - faktor * Ag(k,k:n);
        bg(i) = bg(i) - faktor * bg(k);
    end
end

for i = n:-1:1
    jumlah = 0;
    for j = i+1:n
        jumlah = jumlah + Ag(i,j) * xGauss(j);
    end
    xGauss(i) = (bg(i) - jumlah) / Ag(i,i);
end
waktuGauss = toc;

% METODE GAUSS-JORDAN
Aj = A;
bj = b;

tic;
for k = 1:n
    pivot = Aj(k,k);
    Aj(k,:) = Aj(k,:) / pivot;
    bj(k) = bj(k) / pivot;

    for i = 1:n
        if i ~= k
            faktor = Aj(i,k);
            Aj(i,:) = Aj(i,:) - faktor * Aj(k,:);
            bj(i) = bj(i) - faktor * bj(k);
        end
    end
end

xGJ = bj;
waktuGJ = toc;

% METODE DEKOMPOSISI LU (DOOLITTLE)
L = zeros(n,n);
U = zeros(n,n);
y = zeros(n,1);
xLU = zeros(n,1);

tic;

% Dekomposisi A = L*U
for i = 1:n
    L(i,i) = 1;

    for j = i:n
        jumlah = 0;
        for k = 1:i-1
            jumlah = jumlah + L(i,k) * U(k,j);
        end
        U(i,j) = A(i,j) - jumlah;
    end

    for j = i+1:n
        jumlah = 0;
        for k = 1:i-1
            jumlah = jumlah + L(j,k) * U(k,i);
        end
        L(j,i) = (A(j,i) - jumlah) / U(i,i);
    end
end

% Substitusi maju: L*y = b
for i = 1:n
    jumlah = 0;
    for j = 1:i-1
        jumlah = jumlah + L(i,j) * y(j);
    end
    y(i) = (b(i) - jumlah) / L(i,i);
end

% Substitusi mundur: U*x = y
for i = n:-1:1
    jumlah = 0;
    for j = i+1:n
        jumlah = jumlah + U(i,j) * xLU(j);
    end
    xLU(i) = (y(i) - jumlah) / U(i,i);
end

waktuLU = toc;

% METODE JACOBI
xJ = zeros(n,1);
riwayatJ = zeros(maks_iterasi,n);
galatRelJ = zeros(maks_iterasi,1);
galatRelJVar = zeros(maks_iterasi,n);

tic;
for kJ = 1:maks_iterasi
    xJbaru = zeros(n,1);

    xJbaru(1) = (-1 + xJ(2) + xJ(3)) / 5;
    xJbaru(2) = (2 + xJ(1) + xJ(3) + xJ(4)) / 5;
    xJbaru(3) = (6 + xJ(1) + xJ(2) + xJ(4) + xJ(5)) / 4;
    xJbaru(4) = (2 - xJ(3) + 2*xJ(5)) / 4;
    xJbaru(5) = (-1 - xJ(2) + xJ(3) - xJ(4)) / 4;

    galatRelJVar(kJ,:) = ...
        abs((xJbaru - xJ) ./ xJbaru)';
    galatRelJ(kJ) = max(galatRelJVar(kJ,:));

    riwayatJ(kJ,:) = xJbaru';
    xJ = xJbaru;

    if galatRelJ(kJ) < toleransi
        break;
    end
end
waktuJacobi = toc;
iterJacobi = kJ;

% METODE GAUSS-SEIDEL
xGS = zeros(n,1);
riwayatGS = zeros(maks_iterasi,n);
galatRelGS = zeros(maks_iterasi,1);
galatRelGSVar = zeros(maks_iterasi,n);

tic;
for kGS = 1:maks_iterasi
    xLama = xGS;

    xGS(1) = (-1 + xLama(2) + xLama(3)) / 5;

    xGS(2) = (2 + xGS(1) + ...
              xLama(3) + xLama(4)) / 5;

    xGS(3) = (6 + xGS(1) + xGS(2) + ...
              xLama(4) + xLama(5)) / 4;

    xGS(4) = (2 - xGS(3) + 2*xLama(5)) / 4;

    xGS(5) = (-1 - xGS(2) + ...
              xGS(3) - xGS(4)) / 4;

    galatRelGSVar(kGS,:) = ...
        abs((xGS - xLama) ./ xGS)';
    galatRelGS(kGS) = max(galatRelGSVar(kGS,:));

    riwayatGS(kGS,:) = xGS';

    if galatRelGS(kGS) < toleransi
        break;
    end
end
waktuSeidel = toc;
iterSeidel = kGS;

fprintf('\n========== HASIL METODE GAUSS ==========\n');
fprintf('p = %.10f\n', xGauss(1));
fprintf('q = %.10f\n', xGauss(2));
fprintf('r = %.10f\n', xGauss(3));
fprintf('s = %.10f\n', xGauss(4));
fprintf('t = %.10f\n', xGauss(5));
fprintf('Running time = %.10f detik\n', waktuGauss);

fprintf('\n========== HASIL GAUSS-JORDAN ==========\n');
fprintf('p = %.10f\n', xGJ(1));
fprintf('q = %.10f\n', xGJ(2));
fprintf('r = %.10f\n', xGJ(3));
fprintf('s = %.10f\n', xGJ(4));
fprintf('t = %.10f\n', xGJ(5));
fprintf('Running time = %.10f detik\n', waktuGJ);

fprintf('\n========== HASIL DEKOMPOSISI LU ==========\n');
fprintf('p = %.10f\n', xLU(1));
fprintf('q = %.10f\n', xLU(2));
fprintf('r = %.10f\n', xLU(3));
fprintf('s = %.10f\n', xLU(4));
fprintf('t = %.10f\n', xLU(5));
fprintf('Running time = %.10f detik\n', waktuLU);

fprintf('\n========== HASIL JACOBI ==========\n');
fprintf('Jumlah iterasi = %d\n', iterJacobi);
fprintf('Galat relatif maksimum akhir = %.10f\n', galatRelJ(iterJacobi));
fprintf('Running time = %.10f detik\n', waktuJacobi);
fprintf('Solusi [p q r s t]:\n');
disp(xJ');

fprintf('\nRiwayat iterasi Jacobi [p q r s t]:\n');
disp(riwayatJ(1:iterJacobi,:));

fprintf('Galat relatif Jacobi, iterasi 1 dan 2:\n');
disp(100 * galatRelJVar(1:min(2,iterJacobi),:));

fprintf('\n========== HASIL GAUSS-SEIDEL ==========\n');
fprintf('Jumlah iterasi = %d\n', iterSeidel);
fprintf('Galat relatif maksimum akhir = %.10f\n', galatRelGS(iterSeidel));
fprintf('Running time = %.10f detik\n', waktuSeidel);
fprintf('Solusi [p q r s t]:\n');
disp(xGS');

fprintf('\nRiwayat iterasi Gauss-Seidel [p q r s t]:\n');
disp(riwayatGS(1:iterSeidel,:));

fprintf('Galat relatif Gauss-Seidel, iterasi 1 dan 2:\n');
disp(100 * galatRelGSVar(1:min(2,iterSeidel),:));

% GALAT ABSOLUT METODE ITERATIF(Hampiran) DENGAN METODE GAUSS(Eksak)
galatAbsolutJ = abs(xGauss - xJ);
galatAbsolutGS = abs(xGauss - xGS);

fprintf('\n========== GALAT ABSOLUT ==========\n');
fprintf('Variabel     Jacobi        Gauss-Seidel\n');

nama = {'p','q','r','s','t'};
for i = 1:n
    fprintf('%s       %.10f    %.10f\n', ...
        nama{i}, galatAbsolutJ(i), galatAbsolutGS(i));
end

% RUNNING TIME
fprintf('\n========== RINGKASAN RUNNING TIME ==========\n');
fprintf('Gauss          : %.10f detik\n', waktuGauss);
fprintf('Gauss-Jordan   : %.10f detik\n', waktuGJ);
fprintf('Dekomposisi LU : %.10f detik\n', waktuLU);
fprintf('Jacobi         : %.10f detik\n', waktuJacobi);
fprintf('Gauss-Seidel   : %.10f detik\n', waktuSeidel);

