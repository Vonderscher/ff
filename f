


import matplotlib.pyplot as plt
import numpy as np
Conso_bact=4.e+5 #molecule O2 /s
Avo=6.022e+23 #molecule/mol
R_bact=1.e-6/2. #m
Vol_bact=4*3.14*(R_bact**3.)/3 #m3
phi_bact=0.3
c_o2=7. #mol.m-3
e=500.e-6 #m
D=1.2e-9 #m2.s-1

#Consommation O2 par unité de volume de bactérie en mol.m-3.s-1
cvbact=Conso_bact/(Avo*Vol_bact)
#Consommation O2 par unité de volume de biofilm en mol.m-3.s-1
cvbiof=cvbact*phi_bact
#Constante de vitesse de consommation de l'O2 dans le biofilm s-1
kr=cvbiof/c_o2
Ha=e*np.sqrt(kr/D)
print ('Hatta=',Ha)

jhat=Ha*np.tanh(Ha)
print('jhat =',jhat, '-')
c_i=7
j=jhat*D*c_i/e
print ('j    =',j, 'mol.m-2.s-1')
z=np.linspace(0.,1.,100)
c=np.cosh(Ha*(z-1))/np.cosh(Ha)
plt.plot(z,c, 'ro-')
plt.annotate('', xy=(0.5, 0.5), xytext=(0., 0.5),arrowprops=dict(arrowstyle="->"))
plt.text(0., 0.55, 'mass flux = %.2E mol.m-2.s-1'%j)
plt.xlabel('z/e')
plt.ylabel('c/c_b')
plt.show()

u=1. #m.s-1
d=0.001 #m
ro=1000. #kg.m-3
mu=0.001 #kg.m-1.s-1
Dw=2.e-9 #m2.s-1

re=ro*u*d/mu
sc=mu/(ro*Dw)
print('Reynolds : ', re)
if re<2000 :
    L=1 #m
    sh=1.86*((re*sc*d/L)**0.33)
else :
    sh=0.023*(re**0.8)*(sc**0.33)
print ('Sherwood : ', sh)
k=sh*Dw/d
delta=d/sh
print ('Coefficient de transfert de matière : ', k, 'm/s')
print ('Epaisseur de couche limite massique : ', delta, 'm')

Bi=Dw*e/(D*delta)
print ('Biot :', Bi)
chat=Bi/(Bi+Ha*np.tanh(Ha))
c_b=c_i/chat
print ('Concentration à l interface, ci =', c_i, '    mol.m-3')
j_f=k*(c_b-c_i)
j_b=D*c_i*Ha*np.tanh(Ha)/e
print ('Flux d oxygene consomme      j  =', j_f, 'mol.m-2.s-1')
import matplotlib.patches as mpatches
z=np.linspace(-delta,e,100)
c=c_i*np.cosh(Ha*((z/e)-1))/np.cosh(Ha)
c=np.where(z > 0., c, c_i+(c_i-c_b)*z/delta) 
jz=-(D*c_i/e)*Ha*np.sinh(Ha*((z/e)-1))/np.cosh(Ha)
jz=np.where(z > 0., jz, k*(c_b-c_i)) 
print ('Profil de concentration dans la couche limite puis dans le biofilm')
fig,ax = plt.subplots(1)
plt.plot(z,c, 'ro-')
plt.annotate('', xy=(0, c_b/2), xytext=(-delta/2, c_b/2),arrowprops=dict(arrowstyle="fancy",color='blue'))
plt.text(-delta, 0.05+c_b/2, 'mass flux = %.2E mol.m-2.s-1'%j_f,color='blue')
plt.annotate('', xy=(delta/2,c_b/2), xytext=(0, c_b/2),arrowprops=dict(arrowstyle="fancy",color='green'))
plt.text(0, c_b/2-0.5, 'mass flux = %.2E mol.m-2.s-1'%j_b,color='green')
rectb = mpatches.Rectangle([0, 0], e, c_b, fc="palegreen")
rectf = mpatches.Rectangle([0, 0], -delta, c_b, fc="paleturquoise")
ax.add_patch(rectb)
ax.add_patch(rectf)
plt.xlabel('z')
plt.ylabel('c/c_b')
plt.show()
print ('Densité de flux d oxygene (permet de vérifier la continuité du flux à l interface)')
fig2,ax = plt.subplots(1)
plt.plot(z,jz, 'ro-')
rectb = mpatches.Rectangle([0, 0], e, c_b, fc="palegreen")
rectf = mpatches.Rectangle([0, 0], -delta, c_b, fc="paleturquoise")
ax.add_patch(rectb)
ax.add_patch(rectf)
plt.xlabel('z')
plt.ylabel('j')
plt.show()
