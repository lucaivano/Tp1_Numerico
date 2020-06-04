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
        lambdas[0].append(n + 2)
        lambdas[1].append(lambda_n)
    return lambdas

# ---------------------------------------------------------------

# ------------------------- Orden de Convergencia -------------------------------

def estimarOrdenConvergencia(historiaRaices, nIteraciones) :  
    
    alfa= np.zeros((nIteraciones-1,2))
    
    for n in range(3-1,nIteraciones -1):
        e_n_mas_1 = historiaRaices[n+1][1]-historiaRaices[n][1]
        e_n = historiaRaices[n][1]-historiaRaices[n-1][1]
        e_n_menos_1 = historiaRaices[n-1][1]-historiaRaices[n-2][1]
        
        if((np.abs(e_n_mas_1/e_n) <= np.abs(e_n/e_n_menos_1)) or n<=3): #en el informe se explica el porqué de esta condición
            alfa[n]= n,np.abs(np.log10(np.abs(e_n_mas_1/e_n))/np.log10(np.abs(e_n/e_n_menos_1))) 
        else:
            while(n < nIteraciones-1):
                alfa[n] = n, alfa[n-1][1]
                n=n+1     
            return alfa
           
    return alfa

# ----------------------------------------------------------------

# --------------------Parámetros de la Configuracion -------------------------

tolerancia = 1e-5
#tolerancia = 1e-13
max_it = 100

#Límites del intervalo para bisección
a = 0
b = 2

#Semilla para métodos de Newton Raphhson
x0 = 1.0

#Semillas para método de Secante
sem_sec_0 = 0
sem_sec_1 = 2

#Órdenes de convergencia teóricos para utilizar en el cálculo de lambda
orden_biseccion = 1
orden_secante = 1.618
orden_NR = 2
orden_modif = 2

#---------------------------------------------- Función 1 ------------------------------------------------------------------

raiz_Biseccion,nIteraciones_Biseccion,historiaRaices_Biseccion = biseccion(funcion1,a,b,tolerancia,max_it)
raiz_NR_norm,nIteraciones_NR_norm,historiaRaices_NR_norm = NR_normal(funcion1,x0,der_prim_f1,tolerancia,max_it)
raiz_NR_modif,nIteraciones_NR_modif,historiaRaices_NR_modif = NR_modif(funcion1,x0,der_prim_f1,der_seg_f1,tolerancia,max_it)
raiz_Secante,nIteraciones_Secante,historiaRaices_Secante = Secante(funcion1 ,sem_sec_0, sem_sec_1, tolerancia,max_it)

print("\nBisección")
print_result(historiaRaices_Biseccion, tolerancia)
print("\nNewton Raphson")
print_result(historiaRaices_NR_norm, tolerancia)
print("\nNewton Raphson modificado")
print_result(historiaRaices_NR_modif, tolerancia)
print("\nSecante")
print_result(historiaRaices_Secante, tolerancia)

"""#Descomentar para correr métodos para la función 2
#---------------------------------------------- Función 2 ------------------------------------------------------------------

raiz_Biseccion,nIteraciones_Biseccion,historiaRaices_Biseccion = biseccion(funcion2,a,b,tolerancia,max_it)
raiz_NR_norm,nIteraciones_NR_norm,historiaRaices_NR_norm = NR_normal(funcion2,x0,der_prim_f2,tolerancia,max_it)
raiz_NR_modif,nIteraciones_NR_modif,historiaRaices_NR_modif = NR_modif(funcion2,x0,der_prim_f2,der_seg_f2,tolerancia,max_it)
raiz_Secante,nIteraciones_Secante,historiaRaices_Secante = Secante(funcion2 ,sem_sec_0, sem_sec_1, tolerancia,max_it)

print("\nBisección")
print_result(historiaRaices_Biseccion, tolerancia)
print("\nNewton Raphson")
print_result(historiaRaices_NR_norm, tolerancia)
print("\nNewton Raphson modificado")
print_result(historiaRaices_NR_modif, tolerancia)
print("\nSecante")
print_result(historiaRaices_Secante, tolerancia)
"""

""" #Descomentar para correr métodos para la función 3
#---------------------------------------------- Función 3 ------------------------------------------------------------------

#Utilizar estas semillas para ver la historia de cada método para la función 3 ya que con las originales no convergen los métodos de Newton Rapshon, Newton Raphson 
#modificado y Secante
#x0 = 1.3
#sem_sec_0 = 0.8

raiz_Biseccion,nIteraciones_Biseccion,historiaRaices_Biseccion = biseccion(funcion3,a,b,tolerancia,max_it)

raiz_NR_norm,nIteraciones_NR_norm,historiaRaices_NR_norm = NR_normal(funcion3,x0,der_prim_f3,tolerancia,max_it)

raiz_NR_modif,nIteraciones_NR_modif,historiaRaices_NR_modif = NR_modif(funcion3,x0,der_prim_f3,der_seg_f3,tolerancia,max_it)

raiz_Secante,nIteraciones_Secante,historiaRaices_Secante = Secante(funcion3 ,sem_sec_0, sem_sec_1, tolerancia,max_it)
"""

""" #Descomentar para ver las historias de cada método con las semillas que propusimos con las que convergen los últimos 3 métodos para la función 3
print("\nBisección")
print_result(historiaRaices_Biseccion, tolerancia)
print("\nNewton Raphson")
print_result(historiaRaices_NR_norm, tolerancia)
print("\nNewton Raphson modificado")
print_result(historiaRaices_NR_modif, tolerancia)
print("\nSecante")
print_result(historiaRaices_Secante, tolerancia)
"""


orden_conv_biseccion = estimarOrdenConvergencia(historiaRaices_Biseccion,nIteraciones_Biseccion)
orden_conv_NR = estimarOrdenConvergencia(historiaRaices_NR_norm,nIteraciones_NR_norm)
orden_conv_NRmodif = estimarOrdenConvergencia(historiaRaices_NR_modif,nIteraciones_NR_modif)
orden_conv_secante = estimarOrdenConvergencia(historiaRaices_Secante,nIteraciones_Secante)

lambdas_biseccion = error_asintotico(historiaRaices_Biseccion, orden_biseccion, tolerancia)
lambdas_NR = error_asintotico(historiaRaices_NR_norm, nIteraciones_NR_norm, tolerancia)
lambdas_NRmodif = error_asintotico(historiaRaices_NR_modif, nIteraciones_NR_modif, tolerancia)
lambdas_secante = error_asintotico(historiaRaices_Secante, nIteraciones_Secante, tolerancia)


#Gráfico orden de convergencia
plt.figure()
plt.plot(orden_conv_biseccion[:,0], orden_conv_biseccion[:,1], '--', lw=4, label= 'Bisección')  
plt.plot(orden_conv_secante[:,0], orden_conv_secante[:,1], '-', lw=2, label= 'Secante') 
plt.plot(orden_conv_NR[:,0], orden_conv_NR[:,1], '-', lw=2, label= 'Newton-Raphson') 
plt.plot(orden_conv_NRmodif[:,0], orden_conv_NRmodif[:,1], '-', lw=2, label= 'Newton-Rapson modificado')  

plt.xlabel(r'$n-2$ con $n$:Número de Iteraciones')
plt.ylabel('alfa')
plt.title('orden de convergencia')
plt.legend(loc='best')
plt.grid(True)
plt.show() 

#Gráfico de constantes asintóticas
plt.figure()
plt.plot(lambdas_biseccion[0], lambdas_biseccion[1], '-', lw=2, label= 'Biseccion')
plt.plot(lambdas_secante[0], lambdas_secante[1], '-', lw=2, label= 'secante')
plt.plot(lambdas_NR[0], lambdas_NR[1], '-', lw=2, label= 'NR normal')
plt.plot(lambdas_NRmodif[0], lambdas_NRmodif[1], '-', lw=2, label= 'NR modif')

plt.xlabel(r'$n$ con $n$:Número de Iteraciones')
plt.ylabel('lambda estimado')
plt.title('constante de error asintótica')
plt.legend(loc='best')
plt.grid(True)
plt.show()

"""
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

# ------------------------------------------------------- """
