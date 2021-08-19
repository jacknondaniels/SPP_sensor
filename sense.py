import numpy as np
import cmath as cm
import matplotlib.pyplot as plt

e_0 = 8.85/10**12
c = 299792458

e_1 = 1.7786**2
#au e_2 =  -8.4953 + 1.6j
e_2 = (0.1726 + 3.421j)**2
n = 1.33
e_out = 1
h = 60e-9
h1 = 40e-9
w_1 = 2690930*10**9 #lambda = 700 nm
w_2 = 2975753*10**9 #lambda = 633 nm

def get(i, omega, e_3):
  k_x = cm.sqrt(e_1 * omega**2 / (c**2)) * cm.sin(i)
  k1 = omega/c * cm.sqrt(e_1 - e_1* cm.sin(i))
  Z1 = cm.sqrt(e_1 - e_1* cm.sin(i))/e_1
  glass = np.array([[cm.cos(k1*h), -1j* cm.sin(k1*h)/Z1], [-1j* cm.sin(k1*h)*Z1, cm.cos(k1*h)]])
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
  if ((a.real**2 + a.imag**2) <= 100000000):
    print (k2/e_2 + k3/e_3)
    print (i)
  return (R)


def draw(ax, omega, e_3):
  theta = []
  transmission = []
  rad = 0
  for i in range(1000):
    rad = rad + 0.00157
    theta.append(rad)
    transmission.append(get(rad, omega, e_3))
  ax.plot (theta, transmission, label = str(cm.sqrt(e_3).real))
  min_value = min(transmission)
  min_index = transmission.index(min_value)
  return theta[min_index]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
resonanse = []
refractive = []
for i in range(10):
  resonanse.append(draw (ax1, w_1, n**2))
  refractive.append(n)
  n = n + 0.01
ax1.set_ylabel ("reflection")
ax1.set_xlabel ("Theta, rad")
ax1.legend()
ax2.plot(resonanse, refractive)
ax2.set_ylabel ("index")
ax2.set_xlabel ("resonanse angle, rad")
fig.savefig("hehe6.png")
plt.close()