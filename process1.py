import re
import matplotlib.pyplot as plt

def parse_plot_dat(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    blocks = []
    current_block = None

    for line in lines:
        if line.startswith('j ='):
            if current_block:
                blocks.append(current_block)
            energy = float(re.findall(r'-?\d+\.\d+', line)[0])
            current_block = {'energy': energy, 'x': [], 'y1': [], 'y2': []}
        else:
            values = line.split()
            if len(values) == 3:
                x, y1, y2 = map(float, values)
                current_block['x'].append(x)
                current_block['y1'].append(y1)
                current_block['y2'].append(y2)

    if current_block:
        blocks.append(current_block)

    return blocks

# Ruta al archivo plot.dat
file_path = '/home/joe/Downloads/Mestrado/DIRHB—A relativistic self-consistent mean-field framework for atomic nuclei/Dirhb-package-revised/dirhbs/plotProton.dat'

# Parsear el archivo
#blocks = parse_plot_dat(file_path)

# Imprimir los resultados y graficar
#plt.figure()

# Graficar x vs y1
#plt.figure(figsize=(10, 6))
#for i, block in enumerate(blocks):
#    plt.plot(block['x'], block['y1'], label=f'f 1p 3/2  : { block["energy"]}')
#    plt.plot(block['x'], block['y2'], label=f'g 1p 3/2: {block["energy"]}')
 
#plt.title('Wave function r vs f')
#plt.xlabel('radius')
#plt.ylabel('wave function')
#plt.legend()
#plt.grid(True)
#plt.savefig('snufP.png')
#plt.show()

blocks = parse_plot_dat(file_path)
plt.figure(figsize=(10, 6))
indices = [0, 19, 7, 25, 1, 13, 41, 20, 31, 8]
for i in indices:
    block = blocks[i]
    print(f"Energy: {block['energy']}")
    #plt.plot(block['x'], block['y1'], label=f'Energy: {block["energy"]}')
    plt.plot(block['x'], block['y2'], label=f'Energy: {block["energy"]}')



plt.title('Wave function')
plt.xlabel('radius')
plt.ylabel('wave function')
plt.legend()
plt.grid(True)
plt.savefig(f'dirhbsfPgff.png')
#plt.show()
# Graficar x vs y2
#plt.figure(figsize=(10, 6))
#for j, block in enumerate(blocks):
#    plt.plot(block['x'], block['y2'], label=f'Energy: {block["energy"]}')

#plt.title('Wave function r vs g')
#plt.xlabel('radius')
#plt.ylabel('wave function')
#plt.legend()
#plt.grid(True)
#plt.savefig('snugP.png')
#plt.show()