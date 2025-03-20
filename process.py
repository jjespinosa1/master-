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
file_path = '/home/joe/Downloads/Mestrado/DIRHB—A relativistic self-consistent mean-field framework for atomic nuclei/Dirhb-package-revised/dirhbs/plot1Proton.dat'

# Parsear el archivo
blocks = parse_plot_dat(file_path)

blockn = blocks[13]
print(f"Energy: {blockn['energy']}")
for x, y1, y2 in zip(blockn['x'], blockn['y1'], blockn['y2']):
    #pass  # This line is just to keep the loop structure

# Imprimir los resultados y graficar
#plt.figure(figsize=(10, 6))

    # Graficar 
    plt.plot(blockn['x'], blockn['y1'], label=f'f 1d 3/2  : {blockn["energy"]}')
    plt.plot(blockn['x'], blockn['y2'], label=f'g 1d 7/2: {blockn["energy"]}')
#blocks = parse_plot_dat(file_path)
#plt.figure(figsize=(10, 6))
#indices = [0, 19, 7, 25, 1, 13, 40, 20, 31, 8, 45]
#for i in indices:
#    block = blocks[i]
#    print(f"Energy: {block['energy']}")
    #plt.plot(block['x'], block['y1'], label=f'Energy: {block["energy"]}')
#    plt.plot(block['x'], block['y2'], label=f'Energy: {block["energy"]}')


plt.figure(figsize=(10, 6))
plt.title('Wave function')
plt.xlabel('radius')
plt.ylabel('wave function')
plt.legend()
plt.grid(True)
plt.savefig(f'dirhbsfP-{blockn["energy"]}.png')
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
#plt.savefig('dirhbhgP.png')
#plt.show()