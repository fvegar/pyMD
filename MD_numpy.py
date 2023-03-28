# DINAMICA MOLECULAR EN python, DISCOS DUROS    #####
# Fisica Estadistica, UEx, en Roma, mayo 2015   #####
                                                #####
# VERSION 1.0                                   #####
                                                #####
# en modo terminal, ejecutar con:               #####
                                                #####
# 'python MD.py'                                #####
                                                #####
#####################################################


# importa librerias necesarias: "math", "numpy, "random", "bisect", "operator" y "system"
# usualmente todas ya vienen instaladas en el nucleo python excepto "numpy" 
# hay una version en git/pyMD
import math
import numpy as np
import random
import bisect # libreria de ordenacion de listas (para lista de cols.)
from operator import itemgetter, attrgetter
# import matplotlib # esta es para usar graficos python. a implementar en nuevas versiones
from os import system, remove
import pandas as pd


# radio de las particulas
R=1.
# tamano del sistema
LX = 10*R #112.09982432795857*R
LY = 10*R #112.09982432795857*R
# tamano del sistema menos un radio (para situar las particulas)
LXR = LX*0.5-R
LYR = LY*0.5-R
# fraccion de empaquetamiento
#nu = 0.72
# numero de particulas
#npart = int(math.floor(nu*LX*LY/(math.pi*R*R)))

npart = 10
# numero de pasos temporales (cols.)
nt = 100 * npart

# coef. de restitucion
alfa = 1.0
# parametro de control para evitar solapacion de parts. por error numerico 
# es posible que tol=0 funcione
tol = 1.0e-20
# colisiones/part. debe ser real (no borrar nunca el prefactor 1.0)
# es mas, este parametro no debe modificarse
ncp=1.0*nt/npart
# numero de ncps entre snapshots. si icp=1.0*npart/nt -> se guardan todas las cols. 
# se recomienda icp=npart/nt para etapa de desarrollo de codigo
# en todo caso icp < nt  si se quieren guardar datos; icp > nt  si no se quiere guardar

# iteraciones entre snapshots que sale
utermo = 1
#utermo=int(math.ceil(icp))


#   inicializa listas de velocidades y posiciones 
vx = np.array([0. for i in range(npart)])
vy = np.array([0. for i in range(npart)])

x = np.array([0. for i in range(npart)])
y = np.array([0. for i in range(npart)])

#   inicializa listas temporales de T y a2 
temp = np.array([0. for i in range(nt+1)])
a2 = np.array([0. for i in range(nt+1)])

#inicializa listas relacionadas con las colisiones
listacol = []
listacol_orden = []
ij = []

# inicializa el tiempo
t = 0.
dt = 0.
it = 0

#   inicializa el generador aleatorio. cada vez que se lanza la simulacion usa una semilla aleatoria
# es decir, ejecuciones consecutivas hacen simulaciones estadisticamente diferentes (replicas)
# si no se quiere esta propiedad, escribir: random.seed(1)
random.seed()



#########################
##### FUNCION propaga ###
#########################
#   avanza las particulas con v cte un intervalo de tiempo dt

def propaga(dt):
    global vx, vy, x, y
    x = x + vx * dt
    y = y + vy * dt
        
################################
##### fin de FUNCION propaga ###
################################



def midedist(i,j):
    global x, y
    dx = x[i] - x[j]
    dy = y[i] - y[j]
    dist2 = (dx*dx + dy*dy)- 4*R*R
#    if dist2<0:
#        print (i,j)
#        print ("PELIGRO!!!\n")


#########################
##### FUNCION tcol ######
#########################
#   calcula los tiempos de colision p-p. para un par (i,j)

def tcol(i,j):
    global vx, vy, x, y
    dx = x[i] - x[j]
    dy = y[i] - y[j]
    dvx = vx[i] - vx[j]
    dvy = vy[i] - vy[j]
    drdv = dx*dvx + dy*dvy
    # estructura condicional de colision p-p
    # condicion de acercamiento
    if drdv > 0:
        vct = float('inf')
    else:
        dist2 = (dx*dx + dy*dy) - 4*R*R # distancia instantanea entre dos particulas
        raiz=drdv*drdv - dist2 * (dvx*dvx + dvy*dvy) # condicion de solucion real en la condicion de col.
        if raiz < 0:
            vct = float('inf')
            # si hay sol. real, guarda en dt el tiempo de col.
        else:
            vdt = dist2/(math.sqrt(raiz)-drdv)
            # posicion de la colision. si en realidad la colision ocurriria fuera del sistema, descartala
            xicol = x[i] + vx[i]*vdt
            yicol = y[i] + vy[i]*vdt
            xjcol = x[j] + vx[j]*vdt
            yjcol = y[j] + vy[j]*vdt
            # estructura condicional de col. fuera del sistema
            if math.fabs(xicol) > LXR:
                vdt = float('inf')
            elif math.fabs(xjcol) > LXR:
                vdt = float('inf')
            elif math.fabs(yicol) > LYR:
                vdt = float('inf')
            elif math.fabs(yjcol) > LYR:
                vdt = float('inf')
            else:
                # coloca en la lista de colisiones ordenada de menor a mayor
                # usa un algoritmo rapido 'binary search' para la colocacion
                bisect.insort(listacol,[vdt,[i,j]])


################################
##### fin de FUNCION tcol ######
################################




#########################
##### FUNCION tpcol #####
#########################
#   calcula los tiempos de colision particula-muro. para una particula i.
#   identificadores de pared: -1 (izq.), -2 (inf.), -3 (dcha.), -4 (sup.)

def tpcol(i):
    global vx, vy, x, y
    if vx[i] == 0:
        tx = float('inf')
    elif vx[i] < 0:
        ltx = [-(LXR+x[i])/vx[i],-1]
    elif vx[i] > 0:
        ltx = [(LXR-x[i])/vx[i],-3]
        
    if vy[i] == 0:
        ty = float('inf')
    elif vy[i] < 0:
        lty = [-(LYR+y[i])/vy[i],-2]
    elif vy[i] > 0:
        lty = [(LYR-y[i])/vy[i],-4]

    ltm = sorted( [ltx,lty], key=itemgetter(0) )
    vdt = ltm[0][0]
    im = ltm[0][1]
    bisect.insort( listacol, [vdt,[i,im]] )

################################
##### fin de FUNCION tpcol #####
################################



############################
##### FUNCION pcolisiona  ##
############################
# actualiza velocidad de part. que colisiona con pared

def pcolisiona(ii):
    global vx, vy, x, y   
    if ii[1]==-1 or ii[1]==-3:
        vx[ii[0]] = -vx[ii[0]]
    elif ii[1]==-2 or ii[1]==-4:
        vy[ii[0]] = -vy[ii[0]]

##################################
##### fin de FUNCION pcolisiona ##
##################################




############################
##### FUNCION colisiona  ##
############################
# actualiza velocidades de parts. que han colisionado 

def colisiona(par):
    global vx, vy, x, y
    # la 1a particula la llamamos i y la 2a, j
    i=par[0]
    j=par[1]
    
    dx=x[i]-x[j]
    dy=y[i]-y[j]
    
    # construye sigma_ij unitario
    sigma_norma=math.sqrt(dx*dx+dy*dy)
    sigmax=dx/sigma_norma
    sigmay=dy/sigma_norma
    
    # construye g \cdot sigma (g, vel relativa)
    gsigma=(vx[i]-vx[j])*sigmax+(vy[i]-vy[j])*sigmay
    
    # actualiza vel. de 1a. part.
    vx[i]=vx[i]-0.5*(1+alfa)*gsigma*sigmax
    vy[i]=vy[i]-0.5*(1+alfa)*gsigma*sigmay
    
    # actualiza vel. de 2a. part.
    vx[j]=vx[j]+0.5*(1+alfa)*gsigma*sigmax
    vy[j]=vy[j]+0.5*(1+alfa)*gsigma*sigmay

##################################
##### fin de FUNCION colisiona ###
##################################
    

################################
##### FUNCION escribe_estado  ##
################################

# actualiza velocidades de parts. que han colisionado 
## def write_micr_state(ja):
##     global vx, vy, x, y
##     # formatea el nombre de archivo de posiciones y escribelo en disco

##     print ("####### it: ########", it) # imprime it (n. de cols.) #opcional
##     print ("####### no. archivo: ########", ja) # n. de archivo #opcional
##     inum='{0:05d}'.format(ja)

##     nombre='Datos/xy'+inum+'.dat'
##     xy= pd.DataFrame( np.array([[x[i],y[i]] for i in range(len(x))]) )
##     xy.to_csv(nombre, delimiter='\t',\
##                header =['x','y'] , index=False,float_format='%7.3f')
   
##     # formatea el nombre de archivo de posiciones 
##     nombre='Datos/vxvy'+inum+'.dat'
##     vxvy= pd.DataFrame( np.array([[vx[i],vy[i]] for i in range(len(x))]) )
##     vxvy.to_csv(nombre, delimiter='\t', \
##                 header =['vx','vy'] , index=False, float_format='%7.3f')

def write_micr_state(ja):
    global vx, vy, x, y
    # formatea el nombre de archivo de posiciones y escribelo en disco

    print ("####### it: ########", it) # imprime it (n. de cols.) #opcional
    print ("####### no. archivo: ########", ja) # n. de archivo #opcional
    inum='{0:04d}'.format(ja)

    nombre='Datos/xy'+inum+'.dat'
    with open(nombre,'w') as archivo:
        for i in range(npart):
            archivo.write('{0:10.2f} {1:10.2f}\n'.format( x[i], y[i]))
    archivo.closed

    # formatea el nombre de archivo de posiciones 
    nombre='Datos/vxvy'+inum+'.dat'
    with open(nombre,'w') as archivo:
        for i in range(npart):
            archivo.write('{0:10.2f} {1:10.2f}\n'.format( vx[i], vy[i]))
    archivo.closed


###################################
##### FUNCION initialize_random  ##
###################################
    
def initialize_random():
    global vx, vy, x, y
    # colocacion de las particulas
    x[0]=random.uniform(-LXR, LXR)
    y[0]=random.uniform(-LYR, LYR)

    # condicion de solapamiento
    for i in range(1,npart):
       dr=False
       while dr==False:
           x[i]=random.uniform(-LXR, LXR)
           y[i]=random.uniform(-LYR, LYR)
           # condicion de no solapamiento con las pos. ya generadas
           for j in range(0,i):
               dr=((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])>4*R*R)
               if dr==False:
                  break


    #   velocidades aleatorias para las velocidades, distribucion gaussiana
    for i in range(npart):
        vx[i]=np.random.randn()
        vy[i]=np.random.randn()

        
###################################
##### FUNCION calculate_averages ##
###################################

# measures average fields
def calculate_averages(ja):
    global temp, a2, vx, vy, x, y
    temp[ja]=0.
    a2[ja]=0.
    for i in range(npart):
        vv=vx[i]*vx[i]+vy[i]*vy[i]
        temp[ja]=temp[ja]+vv
        a2[ja]=a2[ja]+vv*vv
    temp[ja]=temp[ja]/npart
    a2[ja]=a2[ja]/(temp[ja]*temp[ja]*npart)
    a2[ja]=(a2[ja]-2.0)*0.5


    
#####################################
##### FUNCION write_averages_evol  ##
#####################################

# wites average fields evolution, in a final file
def write_averages_evol():
    nombre='Datos/temp.dat'
    xy= pd.DataFrame( np.array([[temp[i],a2[i]] for i in range(len(temp))]) )
    xy.to_csv(nombre, sep='\t',\
               header =['T','a2'] , index=False,float_format='%8.5f')


##########################################################
##########################################################
#####                                            #########
#####   PROGRAMA PRINCIPAL, con bucle temporal   #########
#####                                            #########
##########################################################
##########################################################


#### INICIALIZACION. coloca particulas, asigna vels. iniciales y calcula t. de cols. iniciales

print("simulacion MD, alfa= ", alfa)
print("cols/part (total): ", ncp)
#print("cols/part entre snapshots: ", ncp)
print("its. entre snapshots: ", utermo)
print("no. de archivos: ", nt/utermo)



#   genera posiciones aleatorias -no solapantes- para las particulas
#   este algoritmo de colocacion es necesario sustituirlo para densidades altas por otro mejor


initialize_random()

write_micr_state(0)

#### bucle en particulas. Calcula tiempos iniciales de colision 

for i in range(npart-1):
    for j in range(i+1,npart):
        tcol(i,j)   # para todos los pares de particulas (i,j) con j>i
for i in range(npart):
    tpcol(i)    # con la pared

it=0


######  inicia bucle temporal principal  (en cols.)  #######

for it in range(1,nt+1):

    # el tiempo mas corto es la col. que de verdad ocurre
    # guarda como tiempo de col. real (1a. componente de 1er elemento de listacol)
    dt=listacol[0][0]*(1-tol)

    # guardamos las etiquetas del par de particulas que colisionan
    # ( los dos elementos de la 2a. componente del primer elemento de listacol)
    ij=listacol[0][1]
    
    # filtra lista de colisiones, eliminando las cols. que ya no ocurriran
    # es decir, se borran los t. de col. calculados para las particulas que colisionaron
    # y por tanto cambiaron su trayectoria
    
    # elimina antiguas colisiones de 1a particula con pares de mayor indice
    listacol=list(filter(lambda x: x[1][0]!=ij[0] , listacol))
    # elimina antiguas colisiones de 1a particula con pares de menor indice
    listacol=list(filter(lambda x: x[1][1]!=ij[0] , listacol))
    
    if ij[1]>0: # si la segunda particula no es un muro:
        # elimina antiguas colisiones de 2a particula con pares de mayor indice
        listalcol=list(filter(lambda x: x[1][0]!=ij[1] , listacol))
        # elimina antiguas colisiones de 2a particula con pares de menor indice
        listacol=list(filter(lambda x: x[1][1]!=ij[1] , listacol))

    
    t=t+dt # actualiza el tiempo fisico (en escala reducida del sistema)

    # y avanza los tiempos de colision, actualizandolos (porque el t ha avanzado dt)
    limit=range(len(listacol))
    c=[[listacol[i][0]-dt,listacol[i][1]] for i in limit]
    listacol=c
    

    # actualiza primero las posiciones de las parts., 
    # justo hasta la colision que primero ocurre
    propaga(dt)


    # actualiza vels. de part(s). que ha(n) colisionado si la col. es con un muro (pcolisiona)
    # la condicion de colision con un muro es que la "segunda particula" tiene indice negativo
    if ij[1]<0:
        pcolisiona(ij)
    # en caso contrario, la col. es entre dos part. (colisiona)
    else:
        colisiona(ij)
    
    # ahora calculamos los tiempos de col. nuevos para las nuevas trayectorias 
    # de las particulas que colisionarion, 
    # las funciones tcol y tpcol ademas recolocaran ordenadamente esos t en listacol
   
    # primera particula
    i=ij[0]
    # nuevos tiempos de col. de la 1a particula que acaba de colisionar
    tpcol(i)    # con la pared
    
    for j in range(i):
        tcol(j,i)   # para todos los pares de particulas (i,j) con j<i
        
    for j in range(i+1,npart):
        tcol(i,j)   # para todos los pares de particulas (i,j) con j>i
    
    # segunda particula, solo si no es un muro (y por tanto, tiene indice positivo)
    if ij[1]>0:
        i=ij[1]
        # nuevos tiempos de col. de la 2a particula que acaba de colisionar
        tpcol(i)    # con la pared
        
        for j in range(i):
            tcol(j,i)   # para todos los pares de particulas (i,j) con j<i
        
        for j in range(i+1,npart):
            tcol(i,j)   # para todos los pares de particulas (i,j) con j>i
    

##  Escribe pos. y vels. iterativamente

    if (it%utermo==0):
        ia=int(it/utermo)
        write_micr_state(ia)
        # medicion de T y a2 
        calculate_averages(ia)

write_averages_evol()
    
######  FIN bucle temporal principal  (en cols.)  #######


##########################################################
##########################################################
#####                                            #########
#####        FIN de PROGRAMA PRINCIPAL           #########
#####                                            #########
##########################################################
##########################################################

