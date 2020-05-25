# -*- coding: utf-8 -*-
"""
Created on Sat May 23 20:17:31 2020

@author: lu Caivano
"""

import numpy as np
import matplotlib.pyplot as plt
#from scipy import arange

# --------------------- Definicion de Funciones de Uso -------------------
def funcion1 (x):
    return((x**2)-2)

def der_prim_f1(x):
    return (2 * x)

def der_seg_f1(x):
    return 2

def funcion2 (x):
    return(x**5 - 6.6 * x**4 + 5.12 * x**3 + 21.312 * x**2 - 38.016 * x + 17.28 )

def der_prim_f2(x):
    return (5 * x**4 - 26.4 * x**3 + 15.36 * x**2 - 42.624 * x + 38.016 )

def der_seg_f2(x):
    return (20 * x**3 - 79.2 * x**2 - 30.72 * x + 42.624)

def funcion3 (x):
    return((x - 1.5) * np.exp(-4 * (x - 1.5)**2))

def der_prim_f3(x):
    return((-8 * x + 12.0) * (x - 1.5) * np.exp(-4 * (x - 1.5)**2) + np.exp(-4 * (x - 1.5)**2))

def der_seg_f3(x):
    return((-24 * x + (-8 * x + 12.0) * (x - 1.5)**2) * np.exp(-4 * (x - 1.5)**2 ))

# --------------------------------------------------
    
# --------------------- Metodos de Busqueda de Raices -----------------------
error_rt ="No existe una raiz en el intervalo dado" 
error_failure= "Procedimiento terminado sin éxito"
result_rt= "La raiz encontrada"

"""
def Biseccion(f,a,b,tol,nmax):
    
    if np.sign(f(a)*f(b))==0:
        print(error_sqrt)
        return None


    i = 1
    FI=f(a) #Arranco comparando con f(a)
   
    
    while i <= N0:
        p=a+(b-a)/2
        
        FP=f(p)
        
        if (FP==0 or (b-a)/2 <TOL):
            return p,i
        
        else:
                i=i+1
            
        if (np.sign(FI*FP)==1):
            a=p
            FI=FP
            
        else:
              b=p
        
    print(error_failure)
    return None 
  """  


def NR_normal(f,sem,f_deriv_prim,tol,nmax):

    historia = np.zeros((nmax, 2))
    p_ant = sem
    #Antes de la primera iteración guardo con índice 0 la semilla
    historia[0] = (0, sem)

    i = 1
    while i < nmax:

        #Genero la sucesión de Newton Raphson
        p = p_ant - f(p_ant)/f_deriv_prim(p_ant)
        #Guardo en la historia el número de iteración con el nuevo p
        historia[i] = (i, p)

        #Si se llegó a la tolerancia deseada retorno la raíz, su número de iteración y la historia
        if np.abs(p - p_ant) < tol:
            #Acorto la historia hasta la última iteración
            historia = historia[:i+1]
            return p, i, historia

        #Si no, actualizo el p anterior e incremento el índice
        p_ant = p
        i += 1
    #Si se llega a la cantidad máxima de iteraciones sin encontrar la raíz, lo aviso y no retorno nada
    print(error_failure)
    return None

def NR_modif(f,sem,f_deriv_seg,tol,nmax):
    return None

def Secante(f,sem0,sem1,tol,nmax):
    return None

# ----------------------------------------------------

# --------------------Parametros de la Configuracion -------------------------

tolerancia = 1e-5
max_it = 50
#a = 0
#b = 2

#raiz_Biseccion,nIteraciones_Biseccion,historiaRaices_Biseccion = Biseccion(f1,a,b,tolerancia,max_it)


x0 = 1.0
raiz_NR_norm,nIteraciones_NR_norm,historiaRaices_NR_norm = NR_normal(funcion1,x0,der_prim_f1,tolerancia,max_it)

print(historiaRaices_NR_norm)

"""
raiz_NR_modif,nIteraciones_NR_modif,historiaRaices_NR_modif = NR_modif(f1,x0,tolerancia,max_it)

sem_sec_0 = 0
sem_sec_1 = 2
raiz_Secante,nIteraciones_Secante,historiaRaices_Secante = Secante(g,p0,tolerancia,max_it)

# --------------------------------------------------------

# -------------------- Graficos de Funciones ----------------------------           
x = arange(0, 3, 0.01)
plt.figure()
plt.plot(x, [ funcion1 (i) for i in x], '-',lw=2,label='funcion1')
plt.plot(x, [ funcion2 (i) for i in x], '-',lw=2,label='funcion2')
plt.plot(x, [ funcion3 (i) for i in x], '-',lw=2,label='funcion3')
plt.xlabel('x')
plt.ylabel('y')
plt.ylim(-2, 2)
plt.xlim(0,2)
plt.title('Graficos de Funciones en el Intervalo de In')  
plt.legend(loc='best')
plt.grid(True)
plt.show()

# -------------------------------------------------------

# ------------------------- Orden de Convergencia -------------------------------


# -------------------------------------------------------
"""