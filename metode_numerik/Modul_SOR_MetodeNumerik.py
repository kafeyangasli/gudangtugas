from sympy.abc import x, y, z
from sympy import Eq, Poly
from tabulate import tabulate

def SOR_3P(pers_1, pers_2, pers_3, omega, error):
   ## Fungsi Pembantu
   def koefisien(persamaan):
      return (Poly(((persamaan.subs(x, x**3)).subs(y, x**2)).subs(z, x)).all_coeffs())[0:3]

   ## Algoritma Utama
   iterasi = 0                                              # Penunjuk saat ini sudah sampai iterasi keberapa
   hasil = [[0, 0, 0, 0]]                                   # Inisialisasi tabel hasil komputasi
   matriks_a = [                                                           
      koefisien(pers_1.lhs),                                #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
      koefisien(pers_2.lhs),                                # Inisialisasi Matriks A yang merupakan koefisien dari persamaan  #
      koefisien(pers_3.lhs)                                 #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
   ]   
   matriks_x = [x, y, z]                                    # Inisialisasi Matriks X yang berisikan variabel (x, y, z)
   matriks_b = [pers_1.rhs, pers_2.rhs, pers_3.rhs]         # Inisialisasi Matriks B yang merupakan hasil dari persamaan
   diagonal_utama = [matriks_a[i][i] for i in range(3)]     # Elemen dari diagonal utama Matriks A

   persamaan_baru = [                                       # Menginisialisasi persamaan iterasi Gauss-Seidel
      (
         1 / diagonal_utama[i]
      ) * (
         matriks_b[i] 
         - sum(matriks_a[i][j]*matriks_x[j] for j in range(0, i + 1) if j < i) 
         - sum(matriks_a[i][j]*matriks_x[j] for j in range(i, 3) if j > i)
      ) for i in range(3)
   ]

   ### Inisialisasi Iterasi Pertama
   x1 = persamaan_baru[0].subs(y, hasil[0][2]).subs(z, hasil[0][3]).evalf()
   y1 = persamaan_baru[1].subs(x, x1).subs(z, hasil[0][3]).evalf()
   z1 = persamaan_baru[2].subs(x, x1).subs(y, y1).evalf()
   ### Relaksasi Hasil Iterasi Pertama
   x1 = (1 - omega) * hasil[0][1] + omega * x1
   y1 = (1 - omega) * hasil[0][2] + omega * y1
   z1 = (1 - omega) * hasil[0][3] + omega * z1

   iterasi += 1
   hasil.append([iterasi, x1, y1, z1])

   while ( abs(hasil[iterasi][1] - hasil[iterasi-1][1]) > error
        or abs(hasil[iterasi][2] - hasil[iterasi-1][2]) > error
        or abs(hasil[iterasi][3] - hasil[iterasi-1][3]) > error):
      ### Iterasi Gauss-Seidel ke-i
      xn = persamaan_baru[0].subs(y, hasil[iterasi][2]).subs(z, hasil[iterasi][3]).evalf()
      yn = persamaan_baru[1].subs(x, xn).subs(z, hasil[iterasi][3]).evalf()
      zn = persamaan_baru[2].subs(x, xn).subs(y, yn).evalf()
      ### Relaksasi Hasil Iterasi ke-i
      xn = (1 - omega) * hasil[iterasi][1] + omega * xn
      yn = (1 - omega) * hasil[iterasi][2] + omega * yn
      zn = (1 - omega) * hasil[iterasi][3] + omega * zn

      iterasi += 1
      hasil.append([iterasi, xn, yn, zn])
   
   return hasil[-1][1:], tabulate(hasil, headers=["i", "x", "y", "z"], tablefmt='html')