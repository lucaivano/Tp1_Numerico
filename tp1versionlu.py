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
    return(x**5 - 6.6 * x**4 + 5.12 * x**3 + 21.312 * x**2 - 38.016 * x + 17.28 )

def der_prim_f2(x):
    return (5 * x**4 - 26.4 * x**3 + 15.36 * x**2 + 42.624 * x - 38.016 )

def der_seg_f2(x):
    return (20 * x**3 - 79.2 * x**2 + 30.72 * x + 42.624)

def funcion3 (x):
    return((x - 1.5) * np.exp(-4 * (x - 1.5)**2))

def der_prim_f3(x):
    return((-8 * x + 12.0) * (x - 1.5) + 1)* np.exp(-4 * (x - 1.5)**2)

def der_seg_f3(x):
    return((-24 * x + ((8 * x - 12.0)**2) * (x - 1.5)+ 36.0) * np.exp(-4 * (x - 1.5)**2 ))

# --------------------------------------------------

#función para imprimir la tabla de iteraciones indicando su tolerancia 
SEPARADOR_RENGLON = "----------------------------------------------------------------------"  
SEPARADOR_CAMPOS = "\t\t\t" #si queremos exportarlo como CSV 

def print_result(historia,tolerancia):
    
    print(SEPARADOR_RENGLON)
    print('Iteración', "\t\t",'Raíz')
    print(SEPARADOR_RENGLON)

    for elemento in historia: 
        print(elemento[0], SEPARADOR_CAMPOS, elemento[1])
    
    print(SEPARADOR_RENGLON)
    print('con tolerancia = ',tolerancia)
    print(SEPARADOR_RENGLON)

# --------------------- Metodos de Busqueda de Raices -----------------------

ERROR_MSG = "No se encontró una raíz en el intervalo ingresado, seleccione otro por favor."  
NUM_DATOS_TABULADOS = 2  
EXIT_FAILURE = "El procedimiento no convergió." 

def biseccion(f, a, b, tol, n_max):

    if np.sign(f(a)) * np.sign(f(b)) >= 0:
        print(ERROR_MSG)
        return None

    historia = np.zeros((n_max, NUM_DATOS_TABULADOS))
    p_anterior = a
    historia[0] = (0, p_anterior)
    i = 1

    while i <= n_max:

        p = (a + b) / 2
        historia[i] = (i, p)

        if np.abs(p - p_anterior) < tol:
            historia = historia[: i + 1]
            return p, i, historia

        if np.sign(f(a)) * np.sign(f(p)) > 0:
            a = p

        else:
            b = p

        i += 1
        p_anterior = p

    print(EXIT_FAILURE)

    return None

def NR_normal(f,sem,f_deriv_prim,tol,nmax):

    historia = np.zeros((nmax, NUM_DATOS_TABULADOS))
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
    print(EXIT_FAILURE)
    return None

def NR_modif(f,sem,f_deriv_prim,f_deriv_seg,tol,nmax):

    #Defino u(x) y su derivada para aplicarle a esta función el método de Newton Raphson
    def u(x):
        return (f(x)/f_deriv_prim(x))
    def der_u(x):
        return ((f_deriv_prim(x)*f_deriv_prim(x)-f(x)*f_deriv_seg(x))/(f_deriv_prim(x)*f_deriv_prim(x)))
    
    return NR_normal(u, sem, der_u, tol, nmax)


def Secante(f,sem0,sem1,tol,nmax):
    historia = np.zeros((nmax, NUM_DATOS_TABULADOS))
    p_x0 = sem0
    p_x1 = sem1
    #Antes de la primera iteración guardo con índice 0 la semilla
    historia[0] = (0, p_x0)
    historia[1] = (1, p_x1)

    i = 2
    while i < nmax:

        #Genero la sucesión de Secante 
        p = p_x1 - ((f(p_x1) * (p_x1 - p_x0))/(f(p_x1) - f(p_x0)))
        #Guardo en la historia el número de iteración con el nuevo p
        historia[i] = (i, p)

        #Si se llegó a la tolerancia deseada retorno la raíz, su número de iteración y la historia
        if np.abs(p - p_x1) < tol:
            #Acorto la historia hasta la última iteración
            historia = historia[:i+1]
            return p, i, historia

        #Si no, actualizo el p anterior e incremento el índice
        p_x0 = p_x1
        p_x1 = p
        i += 1
    #Si se llega a la cantidad máxima de iteraciones sin encontrar la raíz, lo aviso y no retorno nada
    print(EXIT_FAILURE)
    return None

# ----------------------------------------------------

# ----------------------- Constante De Error Asintótico -------------------------

def error_asintotico(historia, orden, tolerancia):
    diferencias = [abs(historia[n + 1][1] - historia[n][1]) for n in range(len(historia) - 1)]
    lambdas = [[],[]]
    for n in range(len(diferencias) - 1):
        if diferencias[n + 1] < tolerancia:
            break
        lambda_n = diferencias[n + 1] / diferencias[n] ** orden
        lambdas[0].append(n + 3)
        lambdas[1].append(lambda_n)
    return lambdas

# ---------------------------------------------------------------

# ------------------------- Orden de Convergencia -------------------------------

def estimarOrdenConvergencia(historia):
    diferencias = [abs(historia[n + 1][1] - historia[n][1]) for n in range(len(historia) - 1)]
    cocientes = [diferencias[m + 1] / diferencias[m] for m in range(len(diferencias) - 1)]
    ordenes = [[2],[0]]
    for n in range(len(cocientes) - 1):
        if cocientes[n + 1] > cocientes[n]:
            orden_n = ordenes[1][-1] #En lugar de un orden que no tendría sentido usamos el anterior
        else:
            orden_n = np.log(cocientes[n + 1]) / np.log(cocientes[n])
        ordenes[0].append(n + 3)
        ordenes[1].append(orden_n)
    return ordenes

# ----------------------------------------------------------------

# --------------------Parametros de la Configuracion -------------------------

tolerancia = 1e-5
max_it = 500
a = 0
b = 2

raiz_Biseccion,nIteraciones_Biseccion,historiaRaices_Biseccion = biseccion(funcion1,a,b,tolerancia,max_it)


x0 = 1.0
raiz_NR_norm,nIteraciones_NR_norm,historiaRaices_NR_norm = NR_normal(funcion2,x0,der_prim_f2,tolerancia,max_it)

raiz_NR_modif,nIteraciones_NR_modif,historiaRaices_NR_modif = NR_modif(funcion2,x0,der_prim_f2,der_seg_f2,tolerancia,max_it)


sem_sec_0 = 0
sem_sec_1 = 2
raiz_Secante,nIteraciones_Secante,historiaRaices_Secante = Secante(funcion2 ,sem_sec_0, sem_sec_1, tolerancia,max_it)



#print_result(historiaRaices_NR_norm, tolerancia)


'''
# --------------------------------------------------------

# -------------------- Graficos de Funciones ----------------------------           
x = np.arange(0, 3, 0.01)
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
'''
