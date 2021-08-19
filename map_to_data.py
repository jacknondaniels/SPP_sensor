import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
from math import fabs

e_0 = 8.85/10**12
c = 299792458
e_1 = 1.51**2
e_2 = (0.18344 + 3.4332j)**2 #au for 633 nm
e_out = 1
h = 40e-9

w_1 = 2690930*10**9 #lambda = 700 nm
w_2 = 2975753*10**9 #lambda = 633 nm

def readf(f):
    disp = []
    for line in f:   
        disp.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
    return disp

def dispersion (f, l):
  e_2 = 1 + 1j
  Real = []
  Imaginary = []
  Length = []
  for i in range (len(f)):
    if (fabs(l - f[i][0]*1000) < 100):
      Length.append(f[i][0]*1000)
      Real.append(f[i][1])
      Imaginary.append(f[i][2])
  S_real = np.polyfit (Length, Real, 1)
  rl = l*S_real[0] + S_real[1]
  S_imag = np.polyfit (Length, Imaginary, 1)
  img = l*S_imag[0] + S_imag[1]
  return (rl + img*1j)

def get(i, l, e_3, h1, f):
  omega = 2*3.14*c*10**9/l
  k_x = cm.sqrt(e_1 * omega**2 / (c**2)) * cm.sin(i)
  k1 = omega/c * cm.sqrt(e_1 - e_1* cm.sin(i))
  Z1 = cm.sqrt(e_1 - e_1* cm.sin(i))/e_1
  glass = np.array([[cm.cos(k1*h), -1j* cm.sin(k1*h)/Z1], [-1j* cm.sin(k1*h)*Z1, cm.cos(k1*h)]])
  e_2 = dispersion (f, l)
  k2 = omega/c * cm.sqrt(e_2 - e_1* cm.sin(i))
  Z2 = cm.sqrt(e_2 - e_1* cm.sin(i))/e_2
  metal = np.array([[cm.cos(k2*h), -1j* cm.sin(k2*h)/Z2], [-1j* cm.sin(k2*h)*Z2, cm.cos(k2*h)]])
  k3 = omega/c * cm.sqrt(e_3 - e_1* cm.sin(i))
  Z3 = cm.sqrt(e_3 - e_1* cm.sin(i))/e_3
  probe = np.array([[cm.cos(k3*h1), -1j* cm.sin(k3*h1)/Z3], [-1j* cm.sin(k3*h1)*Z3, cm.cos(k3*h1)]])
  k_out = omega/c * cm.sqrt(e_out - e_1* cm.sin(i))
  Z_out = cm.sqrt(e_out - e_1* cm.sin(i))/e_out
  air = np.array([[cm.cos(k_out*h), -1j* cm.sin(k_out*h)/Z_out], [-1j* cm.sin(k_out*h)*Z_out, cm.cos(k_out*h)]])
  out = air.dot(probe.dot(metal.dot(glass)))
  r = ((out[0][0] + Z_out*out[0][1])*Z1 - (out[1][0] + Z_out*out[1][1]))/((out[0][0] + Z_out*out[0][1])*Z1 + (out[1][0] + Z_out*out[1][1]))
  R = r.real**2 + r.imag**2
  a = k2/e_2 + k3/e_3
  #if ((a.real**2 + a.imag**2) <= 100000000):
  #  print (k2/e_2 + k3/e_3)
  #  print (i)
  return (R)


def draw(l, e_3, h1, disp):
  theta = []
  transmission = []
  resonanse_angle = []
  resonanse = []
  rad = 0
  for i in range(1000): # you need to rewrite algorith of finding the plasmonic resonanse minimum
    rad = rad + 0.00157
    if (rad > 0.55) & (rad < 0.90):
      resonanse_angle.append(rad)
      resonanse.append(get(rad, l, e_3, h1, disp))
  min_value = min(resonanse)
  min_index = resonanse.index(min_value)
  #print (resonanse_angle[min_index], omega, h1, e_3)
  return resonanse_angle[min_index]

def find_sense (l, h1, disp):
  resonanse = []
  refractive = []
  n = 1.33 #probe
  for i in range(10):
    resonanse.append(draw (l, n**2, h1, disp))
    refractive.append(n)
    n = n + 0.01
  S =  np.polyfit(refractive, resonanse, 1)
#  print(fabs(S[0])*180/3.14, l, h1)
  return (fabs(S[0])*180/3.14)

f = open ('data.txt', 'w')
g = open ('disp.txt', 'r')
au = readf(g)
g.close()
lda = np.arange (600, 800, 10)
thick = np.arange(40e-9, 60e-9, 1e-9)
z = np.zeros((len(lda), len(thick)))
for i in range (len(lda)):
  for j in range(len(thick)):
    z[i][j] = find_sense(lda[i], thick[j], au)
    f.write (str(z[i][j]) + ' ')
  f.write ('\n')
  print ("was here")
f.close()
