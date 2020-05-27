

error_sqrt ="No existe una raiz en el intervalo dado" 
error_failure= "Procedimiento terminado sin éxito "
result_sqrt= "La raiz encontrada"

def f(x):
    return x**3 + 4*x**2 - 10

print("f(1) = ", f(1))
print("f(2) = ", f(2))


def biseccion(f,a,b,TOL,N0):

    if np.sign(f(a)*f(b))==0:
        print(error_sqrt)
        return None


    historia = []
    
    i = 1
    FI=f(a) #Arranco comparando con f(a)
   
    
    while i <= N0:
        p=a+(b-a)/2
        historia.append((i,p))
        FP=f(p)
        
        if (FP==0 or (b-a)/2 <TOL):
            return p,i,historia
        
        else:
                i=i+1
            
        if (np.sign(FI*FP)==1):
            a=p
            FI=FP
            
        else:
              b=p
        
    print(error_failure)
    return None 

# Paramentros de la configuracion
tolerancia = 1e-4
a = 1
b = 2
nMax = 1000

raiz,nIteraciones,historia = biseccion(f,a,b,tolerancia,nMax)
print("La raiz encontrada es: ", raiz, " usando N = ", nIteraciones, "iteraciones")
 
print("Tolerancia: ", tolerancia)


print("Tabla con la historia:")


print("Iteración", "\t", "Raiz")

for elemento in historia:
    print(elemento[0], "\t \t", elemento[1])

