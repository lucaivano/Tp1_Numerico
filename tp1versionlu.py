# -*- coding: utf-8 -*-
"""
Created on Sat May 23 20:17:31 2020

@author: lu Caivano
"""

import numpy as np
import matplotlib.pyplot as plt

# --------------------- Definicion de Funciones de Uso -------------------
def funcion1 (x):
    return((x**2)-2)

def der_prim_f1(x):
    return (2 * x)

def der_seg_f1(x):
    return 2

def funcion2 (x):
    return(x**5 - 6.6 * x**4 + 5.12 * x**3 + 21.312 * x**2 - 38.016 * x + 17.82 )

def der_prim_f2(x):
    return (5 * x**4 - 26.4 * x**3 + 15.36 * x**2 - 42.624 * x + 38.016 )

def der_seg_f2(x):
    return (20 * x**3 - 79.2 * x**2 - 30.72 * x + 42.624)

def funcion3 (x):
    return((x - 1.5) * np.exp(-4 * (x - 1.5)**2))

def der_prim_f3(x):
    return((-8 * x + 12.0) * (x - 1.5) * np.exp(-4 * (x - 1.5)**2) + np.exp(-4 * (x - 1.5)**2)

def der_seg_f3(x):
    return((-24 * x + (-8 * x + 12.0) * (x - 1.5)**2) * np.exp(-4 * (x - 1.5)**2 )

# --------------------------------------------------
    
# --------------------- Metodos de Busqueda de Raices -----------------------

def Biseccion(f,a,b,tol,nmax):
    


def NR_normal(f,sem,f_deriv_prim,tol,nmax):


def NR_modif(f,sem,f_deriv_seg,tol,nmax):


def Secante(f,sem0,sem1,tol,nmax):


# ----------------------------------------------------

# --------------------Parametros de la Configuracion -------------------------
tolerancia = 1e-5
a = 0
b = 2

raiz_Biseccion,nIteraciones_Biseccion,historiaRaices_Biseccion = Biseccion(f1,a,b,tolerancia,nMax)

x0 = 1.0
raiz_NR_norm,nIteraciones_NR_norm,historiaRaices_NR_norm = NR_norm(f1,x0,tolerancia,nMax)

raiz_NR_modif,nIteraciones_NR_modif,historiaRaices_NR_modif = NR_modif(f1,x0,tolerancia,nMax)

sem_sec_0 = 0
sem_sec_1 = 2
raiz_Secante,nIteraciones_Secante,historiaRaices_Secante = Secante(g,p0,tolerancia,nMax)

# --------------------------------------------------------

# -------------------- Graficos de Funciones ----------------------------           
x = range(0, 3)
plt.figure()
plt.plot(x, [ funcion1 (i) for i in x], '-',lw=2,label='funcion1')
plt.plot(x, [ funcion2 (i) for i in x], '-',lw=2,label='funcion2')
plt.plot(x, [ funcion3 (i) for i in x], '-',lw=2,label='funcion3')
plt.xlabel('x')
plt.ylabel('y')
plt.ylim(-3, 3)
plt.title('Graficos de Funciones en el Intervalo de In')  
plt.legend(loc='best')
plt.grid(True)
plt.show()

# -------------------------------------------------------

# ------------------------- Orden de Convergencia -------------------------------


# -------------------------------------------------------