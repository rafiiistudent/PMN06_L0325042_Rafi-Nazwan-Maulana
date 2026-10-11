# Penyelesaian Sistem Persamaan Linear (SPL) Menggunakan MATLAB

Program MATLAB ini dibuat untuk menyelesaikan Sistem Persamaan Linear (SPL) dengan lima variabel ($p, q, r, s, t$) menggunakan berbagai metode numerik, baik metode langsung (*direct methods*) maupun metode iteratif (*iterative methods*).

## Metode yang Digunakan

Program ini membandingkan 5 metode berbeda untuk menyelesaikan matriks SPL yang sama:
1. **Metode Gauss** (Eliminasi Gauss)
2. **Metode Gauss-Jordan**
3. **Dekomposisi LU (Metode Doolittle)**
4. **Metode Jacobi** (Iteratif dengan batas toleransi dan maksimal iterasi)
5. **Metode Gauss-Seidel** (Iteratif dengan batas toleransi dan maksimal iterasi)

---

## Fitur & Output Program

Program secara otomatis menghitung dan menampilkan:
* Nilai penyelesaian variabel ($p, q, r, s, t$) untuk masing-masing metode.
* Waktu eksekusi (*running time*) menggunakan fungsi `tic` dan `toc` untuk membandingkan kecepatan tiap metode.
* Riwayat iterasi serta galat relatif maksimum untuk metode iteratif (Jacobi & Gauss-Seidel).
* Perbandingan galat absolut antara metode iteratif terhadap metode Gauss (sebagai solusi eksak).
* Ringkasan performa waktu eksekusi dari kelima metode di akhir program.

---

## Cara Menjalankan Program

1. Pastikan telah menginstal software **MATLAB** di perangkat anda.
2. Unduh atau *clone* repositori ini, lalu buka file `SPL.m`.
3. Jalankan program dengan mengeklik tombol **Run** di MATLAB atau ketik perintah berikut di *Command Window*:
   ```matlab
   SPL
