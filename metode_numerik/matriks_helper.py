def kali(a, b):
   kMA = len(a)
   bMA, bMB = len(a[0]), len(b[0])
   
   r = []
   for i in range(kMA):
      phi = []
      for j in range(bMB):
         ph = 0
         for k in range(bMA):
            ph += a[i][k] * b[k][j]
         phi.append(ph)
      r.append(phi)

   return r

def norm(v):
   return sum([i*i for i in v]) ** 0.5

def id(n):
   return [[0 if i != j else 1 for i in range(n)] for j in range(n)]

def tr(mtx):
   mtx_m = len(mtx)
   mtx_n = len(mtx[0])
   h = []
   for i in range(mtx_n):
      ph = []
      for j in range(mtx_m):
         ph.append(mtx[j][i])
      h.append(ph)
   return h

def roundmtx(mtx):
   m = len(mtx)
   n = len(mtx[0])
   return [[round(mtx[j][i], 5) for i in range(n)] for j in range(m)]

def sign(a):
   return 1 if a >= 0 else -1

def det(mtx):
   m = len(mtx)
   n = len(mtx[0])
   
   if (m != n):
      raise ValueError("Matriks harus NxN!")
   
   if m == 1:
      return mtx[0][0]
   elif m == 2:
      return mtx[0][0] * mtx[1][1] - mtx[0][1] * mtx[1][0]
   elif m == 3:
      a = (mtx[0][0] * mtx[1][1] * mtx[2][2]) + (mtx[0][1] * mtx[1][2] * mtx[2][0]) + (mtx[0][2] * mtx[1][0] * mtx[2][1])
      b = (mtx[0][1] * mtx[1][0] * mtx[2][2]) + (mtx[0][0] * mtx[1][2] * mtx[2][1]) + (mtx[0][2] * mtx[1][1] * mtx[2][0])
      return a - b
   else:
      rv = 0
      for i in range(m):
         submtx = [k[:i] + k[i+1:] for k in mtx[1:]]
         cf = (-1 * i) * mtx[0][i] * det(submtx)
         rv += cf
         
      return rv

def print_matriks(m):
   for baris in m:
      for kolom in baris:
         print(kolom, end=" ")