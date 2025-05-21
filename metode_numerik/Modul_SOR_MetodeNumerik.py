from sympy.abc import x, y, z
from sympy import Eq, Poly
from tabulate import tabulate

class HasilSOR_3P():
   def __init__(self, tabel_hasil, list_persamaan, iterasi, omega, error):
      hasil = tabel_hasil [-1][1:]
      self.relax = omega
      self.error = error
      self.persamaan = list_persamaan
      self.x = hasil[0]
      self.y = hasil[1]
      self.z = hasil[2]
      self.i = iterasi
      self.tabit = tabel_hasil
      self.def_tabit = tabulate(tabel_hasil, headers=["i", "x", "y", "z"])
      self.html_tabit = tabulate(tabel_hasil, headers=["i", "x", "y", "z"], tablefmt="html")

   def __str__(self):
      p1 = self.persamaan[0]
      p2 = self.persamaan[1]
      p3 = self.persamaan[2]
      return f"""HASIL KOMPUTASI SUCCESSIVE OVER-RELAXATION\n{"-" * 10}\n
PERSAMAAN 1\t\t: {str(p1.lhs) + " = " + str(p1.rhs)}
PERSAMAAN 2\t\t: {str(p2.lhs) + " = " + str(p2.rhs)}
PERSAMAAN 3\t\t: {str(p3.lhs) + " = " + str(p3.rhs)}
MENGGUNAKAN GALAT\t: {self.error}
DAN FAKTOR RELAKSASI\t: {self.relax}
{"-" * 10}\n
BANYAKNYA ITERASI\t: {self.i}
NILAI HAMPIRAN X\t: {self.x}
NILAI HAMPIRAN Y\t: {self.y}
NILAI HAMPIRAN Z\t: {self.z}
            """

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
   
   return HasilSOR_3P(hasil, [pers_1, pers_2, pers_3], iterasi, omega, error)

## Contoh Pengaplikasian
if __name__ == "__main__":
   ## SPL SUDAH HARUS TERURUT BERDASARKAN DIAGONAL UTAMA TERBESARNYA!
   ## JIKA TIDAK, MAKA KEMUNGKINAN BESAR AKAN DIVERGEN

   pers_1 = Eq(10*x + 5*y + 2*z, 75)
   pers_2 = Eq(x + 7*y + 5*z, 101)
   pers_3 = Eq(6*x + 3*y + 9*z, 123)

   hasil = SOR_3P(pers_1, pers_2 ,pers_3, [0, 0, 0], 1, 0.01)

   print(hasil)  # Ambil Nilai Hasil Terakhirnya saja
   hasil.html_tabit  # Menampilkan tabel dalam bentuk html agar bisa terlihat jelas dalam IPYNB
   print(hasil.def_tabit) # Atau dalam bentuk print saja