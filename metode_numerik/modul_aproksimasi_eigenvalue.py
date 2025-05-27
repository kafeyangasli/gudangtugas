from numpy import linalg, matrix
from sympy.abc import lamda
from sympy import solve, sqrt

epsilon1 = 10 ** (-2)      # Toleransi galat

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
   
def PowerMethod(A, x0):
   m = len(A)
   n = max(len(A[i]) for i in range(m))
   result = [[0, x0, 0]]
   iterasi = 0

   if (m != n) or (len(x0) != m):
      raise ValueError("Matriks A harus berbentuk nxn dan inisalisasi x0 harus memiliki baris n!")

   while (abs(result[iterasi][2] - result[iterasi-1][2]) > epsilon1 or iterasi == 0):
      yn = tr(kali(A, tr([result[iterasi][1]])))[0]
      norm_yn = max(yn, key=lambda x: abs(x))
      x1 = [i / norm_yn for i in yn]
      iterasi += 1
      result.append([iterasi, x1, norm_yn])

   return result

def InversePowerMethod(A, x0):
   m = len(A)
   n = max(len(A[i]) for i in range(m))
   result = [[0, x0, 0]]
   iterasi = 0

   if (m != n) or (len(x0) != m):
      raise ValueError("Matriks A harus berbentuk nxn dan inisalisasi x0 harus memiliki baris n!")

   q = 0
   for i in range(m):
      for j in range(m):
         q += x0[i] * A[i][j] * x0[j]
   q /= sum(i*i for i in x0)
   qi = id(n); qi = [[qi[i][j] * q for j in range(m)] for i in range(m)]
   aqi = [[A[i][j] - qi[i][j] for j in range(m)] for i in range(m)]

   while (abs(result[iterasi][2] - result[iterasi-1][2]) > epsilon1 or iterasi == 0):
      yn = list(linalg.solve(aqi, result[iterasi][1]))
      norm_yn = max(yn, key=lambda x: abs(x))
      x1 = [i / norm_yn for i in yn]
      iterasi += 1
      result.append([iterasi, [float(i) for i in x1], (1 / norm_yn) + q])
   
   return result

def HouseholderTD(M):
   if (M != tr(M)):
      raise ValueError("Matriks M harus simetris!")

   n = len(M)
   prev_m = M[:]
   nm = list()

   for k in range(0, n - 2):
      alph = (-1 * sign(prev_m[k + 1][k])) * (sum((prev_m[j][k] * prev_m[j][k]) for j in range(k + 1, n)) ** 0.5)
      r = ((0.5 * (alph * alph)) - (0.5 * alph * (prev_m[k + 1][k]))) ** 0.5
      
      if r == 0:
         continue      

      w = [[0 for _ in range(k + 1)] + [(prev_m[k + 1][k] - alph) / (2*r)] + [(prev_m[j][k]) / (2*r) for j in range(k + 2, n)]]
      wt = tr(w)

      wwt = kali(wt, w)
      wwt = [[2 * wwt[i][j] for j in range(n)] for i in range(n)]
      identity = id(n)

      p = [[identity[i][j] - wwt[i][j] for j in range(n)] for i in range(n)]
      
      nm = roundmtx(kali(kali(p, prev_m), p))
      prev_m = nm

   return prev_m

def QR_DECOMP(M):
   m = len(M)
   n = len(M[0])
   q = id(m)
   r = M[:]

   for k in range(min(m, n)):
      a = [r[j][k] for j in range(k, m)]
      norm_a = norm(a)
      sgn_a = sign(a[0])
      e_i = [1] + [0 for _ in range(len(a) - 1)]
      vr = [(sgn_a * norm_a * e_i[j]) for j in range(len(a))] 
      v = [[a[j] + vr[j]] for j in range(len(a))]
      vt = tr(v)
      vvt = kali(v, vt)
      vvt = [[((2 / (norm(vt[0]) ** 2)) * vvt[i][j] if norm(vt[0]) > 0 else 0) for j in range(len(a))] for i in range(len(a))]

      h = id(m)
      for i in range(k, m):
         for j in range(k, m):
            h[i][j] -= ((2 / (norm(vt[0]) ** 2)) if norm(vt[0]) > 0 else 0) * v[i - k][0] * vt[0][j - k]
      h = roundmtx(h)
      
      q = roundmtx(kali(q, tr(h)))
      r = roundmtx(kali(h, r)) 

   return q, r

def QR_EV(M):
   a = M[:]
   for _ in range(10):
      q, r = QR_DECOMP(a)
      a = kali(r, q)
   return a

def SVD(M):
   mt = tr(M)
   mtm = kali(mt, M)
   new_mn = len(mtm)
   identity = id(new_mn)
   identity = [[lamda * identity[i][j] for j in range(new_mn)] for i in range(new_mn)]
   ev_mtx = [[mtm[i][j] - identity[i][j] for j in range(new_mn)] for i in range(new_mn)]
   lamda_eq = det(ev_mtx)
   lambdas = solve(lamda_eq); lambdas = sorted([(lambdas[i].evalf()) for i in range(len(lambdas))], reverse=True)
   diag_mtx = [[lambdas[i] if i == j and i < len(lambdas) else 0 for j in range(len(M[0]) - 1)] for i in range(len(M) - 1)]

   return [sqrt(lambdas[i]) for i in range(len(lambdas))], diag_mtx

