from sympy.abc import x, y, z
from sympy import Eq, Poly
from tabulate import tabulate

def SOR_3P(p1, p2, p3, omega, error):
   ## Fungsi Pembantu
   def koefisien(persamaan):
      return (Poly(((persamaan.subs(x, x**3)).subs(y, x**2)).subs(z, x)).all_coeffs())[0:3]

   ## Algoritma Utama
   ite = 0
   hasil = [[0, 0, 0, 0]]
   matriks_a = [koefisien(p1.lhs), koefisien(p2.lhs), koefisien(p3.lhs)]
   matriks_x = [x, y, z]
   matriks_b = [p1.rhs, p2.rhs, p3.rhs]
   diagonal_utama = [matriks_a[i][i] for i in range(3)]

   persamaan_baru = [
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
   ### Relaksasi
   x1 = (1 - omega) * hasil[0][1] + omega * x1
   y1 = (1 - omega) * hasil[0][2] + omega * y1
   z1 = (1 - omega) * hasil[0][3] + omega * z1

   ite += 1
   hasil.append([ite, x1, y1, z1])

   while ( abs(hasil[ite][1] - hasil[ite-1][1]) > error
        or abs(hasil[ite][2] - hasil[ite-1][2]) > error
        or abs(hasil[ite][3] - hasil[ite-1][3]) > error):
      ### Gauss-Seidel
      xn = persamaan_baru[0].subs(y, hasil[ite][2]).subs(z, hasil[ite][3]).evalf()
      yn = persamaan_baru[1].subs(x, xn).subs(z, hasil[ite][3]).evalf()
      zn = persamaan_baru[2].subs(x, xn).subs(y, yn).evalf()
      ### Relaksasi
      xn = (1 - omega) * hasil[ite][1] + omega * xn
      yn = (1 - omega) * hasil[ite][2] + omega * yn
      zn = (1 - omega) * hasil[ite][3] + omega * zn

      ite += 1
      hasil.append([ite, xn, yn, zn])
   
   return hasil[-1][1:], tabulate(hasil, headers=["i", "x", "y", "z"], tablefmt='html')