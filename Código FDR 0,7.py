import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import tarfile

# === Parámetros del sistema ===
L = 46.606726129198194
N = 484
diametro = 1.0
radio = diametro / 2
phi = 0.7
area_total = L ** 2
densidad = N / area_total

dr = 0.1  # resolución radial

# === Paso 1: Extraer el archivo .tar.gz si aún no está extraído ===
ruta_tar = r"C:\Users\pablo\OneDrive\Escritorio\Universidad Pablo\4 CARRERA\Segundo cuatrimestre\AFES\Trabajo 1\nu07.tar.gz"
ruta_destino = r"C:\Users\pablo\OneDrive\Escritorio\Universidad Pablo\4 CARRERA\Segundo cuatrimestre\AFES\Trabajo 1\nu07_extraido"

if not os.path.exists(ruta_destino):
    os.makedirs(ruta_destino)
    print("Extrayendo archivo...")
    with tarfile.open(ruta_tar, "r:gz") as tar:
        tar.extractall(path=ruta_destino)
    print("Extracción completa.")

# === Paso 2: Obtener los últimos 500 archivos .dat ordenados ===
archivos = sorted(glob.glob(os.path.join(ruta_destino, "*.dat")))
ultimos_500 = archivos[-500:]
print(f"Archivos seleccionados: {len(ultimos_500)}")

# === Función para leer posiciones ===
def leer_posiciones(archivo):
    return np.loadtxt(archivo)

# === Cálculo de g(r) ===
def calcular_gr(posiciones, box_size, dr=0.1, r_max=None):
    N = len(posiciones)
    if r_max is None:
        r_max = box_size / 2

    nbins = int(r_max / dr)
    hist = np.zeros(nbins)

    for i in range(N):
        for j in range(i + 1, N):
            dx = posiciones[i][0] - posiciones[j][0]
            dy = posiciones[i][1] - posiciones[j][1]

            # condiciones de borde periódicas
            dx -= box_size * np.round(dx / box_size)
            dy -= box_size * np.round(dy / box_size)

            r = np.sqrt(dx**2 + dy**2)

            # Si la distancia r es mayor que r_max, se ignora esta pareja
            if r >= r_max:
                continue  # Ignoramos los pares que se salen de r_max

            # Cálculo del índice del bin
            bin_index = int(r / dr)
            
            # Si el índice se pasa, lo ajustamos a nbins-1 (último bin)
            if bin_index >= nbins:
                bin_index = nbins - 1

            hist[bin_index] += 2  # i-j y j-i

    r_vals = (np.arange(nbins) + 0.5) * dr
    area_anillos = 2 * np.pi * r_vals * dr
    norm = N * densidad * area_anillos
    g_r = hist / norm
    return r_vals, g_r


# === Paso 3: Loop para promedio de g(r) ===
gr_acumulado = None
contador = 0

for archivo in ultimos_500:
    posiciones = leer_posiciones(archivo)
    r_vals, g_r = calcular_gr(posiciones, L, dr)

    if gr_acumulado is None:
        gr_acumulado = g_r
    else:
        gr_acumulado += g_r
    contador += 1

g_r_prom = gr_acumulado / contador

# === Paso 4: Graficar ===
plt.figure(figsize=(8, 5))
plt.plot(r_vals, g_r_prom, label="g(r)")
plt.axvline(diametro, color='r', linestyle='--', label='diámetro = 1')
plt.xlabel("r")
plt.ylabel("g(r)")
plt.title("Función de distribución radial2.0 (2D) - últimos 500 archivos .dat")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
