import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Leer el archivo Novo.txt
file_path = '/home/joe/Downloads/Mestrado/DIRHB—A relativistic self-consistent mean-field framework for atomic nuclei/Dirhb-package-revised/dirhbs/Novo1.txt'
#file_path1 = '/home/joe/Downloads/Mestrado/DIRHB—A relativistic self-consistent mean-field framework for atomic nuclei/Dirhb-package-revised/dirhbs/Novo1.txt'
# Leer el archivo en un DataFrame de pandas
df = pd.read_csv(file_path, delim_whitespace=True)

# Extraer las tres primeras columnas
A0 = df['A0']
A1 = df['A1']
A2 = df['A2']

# Crear una figura y un eje 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Graficar los datos
ax.plot(N, Z, A0, c='r')

# Etiquetas de los ejes
ax.set_xlabel('N')
ax.set_ylabel('Z')
ax.set_zlabel('A0')

# Título del gráfico
ax.set_title('Gráfico 3D de N, Z y A0')

# Mostrar el gráfico
plt.show()