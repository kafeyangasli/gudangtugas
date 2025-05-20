from sympy.abc import x, y, z
from sympy import Eq, Poly
from tabulate import tabulate

def SOR_3P(pers_1, pers_2, pers_3, estimate, omega, error):
   ## Fungsi Pembantu
   def koefisien(persamaan):
      return (Poly(((persamaan.subs(x, x**3)).subs(y, x**2)).subs(z, x)).all_coeffs())[0:3]

   ## Algoritma Utama
   iterasi = 0                                              # Penunjuk saat ini sudah sampai iterasi keberapa
   hasil = [[iterasi] + estimate]                           # Inisialisasi tabel hasil komputasi
   matriks_a = [                                                           
      koefisien(pers_1.lhs),                                #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
      koefisien(pers_2.lhs),                                #  Inisialisasi Matriks A yang merupakan koefisien dari persamaan #
      koefisien(pers_3.lhs)                                 #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #
   ]   
   matriks_x = [x, y, z]                                    # Inisialisasi Matriks X yang berisikan variabel (x, y, z)
   matriks_b = [pers_1.rhs, pers_2.rhs, pers_3.rhs]         # Inisialisasi Matriks B yang merupakan hasil dari persamaan

   persamaan_baru = [                                       # Menginisialisasi persamaan iterasi Gauss-Seidel
      (
         1 / matriks_a[i][i]
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
   ### Masukkan hasil iterasi ke tabel sebagai dokumentasi
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
      ### Masukkan hasil iterasi ke tabel sebagai dokumentasi
      iterasi += 1
      hasil.append([iterasi, xn, yn, zn])
   
   return hasil[-1][1:], tabulate(hasil, headers=["i", "x", "y", "z"], tablefmt='html')

## Contoh Pengaplikasian
if __name__ == "__main__":
   ## SPL SUDAH HARUS TERURUT BERDASARKAN DIAGONAL UTAMA TERBESARNYA!
   ## JIKA TIDAK, MAKA KEMUNGKINAN BESAR AKAN DIVERGEN

   pers_1 = Eq(10*x + 5*y + 2*z, 75)
   pers_2 = Eq(x + 7*y + 5*z, 101)
   pers_3 = Eq(6*x + 3*y + 9*z, 123)

   h, t = SOR_3P(pers_1, pers_2 ,pers_3, [0, 0, 0], 1, 0.01)

   h  # Ambil Nilai Hasil Terakhirnya saja
   t  # Menampilkan tabel dalam bentuk html agar bisa terlihat jelas dalam IPYNB