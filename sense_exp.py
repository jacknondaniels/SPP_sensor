import numpy as np
import cmath as cm
import matplotlib.pyplot as plt

e_0 = 8.85/10**12
c = 299792458

e_1 = 1.51**2
e_2 = (0.18344 + 3.4332j)**2 #au for 633 nm
e_out = 1
h = 40e-9
h1 = 20e-9
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
  resonanse_angle = []
  resonanse = []
  rad = 0
  for i in range(1000):
    rad = rad + 0.00157
    theta.append(rad)
    transmission.append(get(rad, omega, e_3))
    if (rad > 0.55) & (rad < 0.90):
      resonanse_angle.append(rad)
      resonanse.append(get(rad, omega, e_3))
  ax.plot (theta, transmission, label = str(cm.sqrt(e_3).real))
  min_value = min(resonanse)
  min_index = resonanse.index(min_value)
  return resonanse_angle[min_index]

def find_sense (h1):
  fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
  resonanse = []
  refractive = []
  n = 1.33 #probe
  for i in range(20):
    resonanse.append(draw (ax1, w_2, n**2))
    refractive.append(n)
    n = n + 0.005
  ax1.set_ylabel ("reflection")
  ax1.set_xlabel ("Theta, rad")
  ax1.legend()
  ax2.plot(refractive, resonanse)
  ax2.set_xlabel ("index")
  ax2.set_ylabel ("resonanse angle, rad")
  fig.savefig("sensetivity"+str(h1)+".png")
  plt.close()
  S =  np.polyfit(refractive, resonanse, 1)
  return (-S[0]*180/3.14)

sensetivity = []
thickness = []
for i in range(10):
  sensetivity.append(find_sense(h1))
  thickness.append(h1)
  h1 = h1 + 4e-9
fig, ax = plt.subplots(1, figsize=(20, 8))
ax.scatter(thickness, sensetivity)
ax.set_ylabel ("sensetivity, deg/RIU")
ax.set_xlabel ("thickness, m")
ax.set_xlim(20e-9, 60e-9)
fig.savefig("thickness of probe layer - sensetivity.png")
plt.close()

