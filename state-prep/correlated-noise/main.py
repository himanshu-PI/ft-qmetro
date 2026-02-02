from funcs import *

def pauli2(p2, x, y):
    str = f'''PAULI_CHANNEL_2({p2/3},0,0, {p2/3},{p2/3},0,0, 0,0,0,0, 0,0,0,0) {x} {y}'''

    return str


def get_ghz_circuit(d,
                     p,
                     q,
                     p2,
                     p_prep, 
                     rounds):

    script = ''
    num_stabs = d - 1

    data_qubits = [i for i in range(d)]

    coords = []
    for i in range(d):
        coords.append(f'QUBIT_COORDS({i},{0},0) {data_qubits[i]}')

    script += '\n'.join(coords) + '\n'


    measure_qubits = [d + i for i in range(num_stabs)]


    coords = []

    for i in range(num_stabs):
        coords.append(f'QUBIT_COORDS({i + 0.5},{0.5},0) {measure_qubits[i]}')

    script += '\n'.join(coords) + '\n'


    script += f"R {' '.join(map(str, data_qubits))}\n"

    script += f"R {' '.join(map(str, measure_qubits))}\n"

    script += f"X_ERROR({p_prep}) {' '.join(map(str, measure_qubits))}\n"



    for i in range(num_stabs):

        script += f'CNOT {data_qubits[i]} {measure_qubits[i]}\n'
        script += pauli2(p2, data_qubits[i], measure_qubits[i]) + '\n'

    script += f"X_ERROR({p}) {data_qubits[-1]}\n"

    for i in range(num_stabs):

        script += f'CNOT {data_qubits[i + 1]} {measure_qubits[i]}\n'
        script += pauli2(p2, data_qubits[i + 1], measure_qubits[i]) + '\n'

    script += f"X_ERROR({p}) {data_qubits[0]}\n"


    script += 'TICK\n'
    script += f"X_ERROR({q}) {' '.join(map(str, measure_qubits))}\n"
    script += f"X_ERROR({p}) {' '.join(map(str, data_qubits))}\n"
    script += f"MR {' '.join(map(str, measure_qubits))}\n"

    round_template = '\t'

    round_template += f"X_ERROR({p_prep}) {' '.join(map(str, measure_qubits))}\n"
    round_template += f"X_ERROR({p}) {' '.join(map(str, data_qubits))}\n\t"



    for i in range(num_stabs):

        round_template += f'CNOT {data_qubits[i]} {measure_qubits[i]}\n\t'
        round_template += pauli2(p2, data_qubits[i], measure_qubits[i]) + '\n\t'

    round_template += f"X_ERROR({p}) {data_qubits[-1]}\n\t"

    for i in range(num_stabs):

        round_template += f'CNOT {data_qubits[i + 1]} {measure_qubits[i]}\n\t'
        round_template += pauli2(p2, data_qubits[i + 1], measure_qubits[i]) + '\n\t'

    round_template += f"X_ERROR({p}) {data_qubits[0]}\n\t"



    round_template += 'TICK\n\t'
    round_template += f"X_ERROR({q}) {' '.join(map(str, measure_qubits))}\n\t"
    round_template += f"MR {' '.join(map(str, measure_qubits))}\n\t"


    for i in range(num_stabs):
        x = -1 - i
        round_template += f"DETECTOR({i + 0.5}, {0}, 0) rec[{x}] rec[{x - num_stabs}]\n\t"

    round_template += 'SHIFT_COORDS(0, 1, 0)\n\t'


    script += f"REPEAT {rounds} {{\n{round_template}}}\n"


    # -----------------------------------

    script += f"X_ERROR({p_prep}) {' '.join(map(str, measure_qubits))}\n"
    script += f"X_ERROR({p}) {' '.join(map(str, data_qubits))}\n"


    for i in range(num_stabs):

        script += f'CNOT {data_qubits[i]} {measure_qubits[i]}\n'
        script += pauli2(p2, data_qubits[i], measure_qubits[i]) + '\n'

    script += f"X_ERROR({p}) {data_qubits[-1]}\n"

    for i in range(num_stabs):

        script += f'CNOT {data_qubits[i + 1]} {measure_qubits[i]}\n'
        script += pauli2(p2, data_qubits[i + 1], measure_qubits[i]) + '\n'

    script += f"X_ERROR({p}) {data_qubits[0]}\n"

    script += 'TICK\n'
    script += f"X_ERROR({q}) {' '.join(map(str, measure_qubits))}\n"
    script += f"MR {' '.join(map(str, measure_qubits))}\n"


    for i in range(num_stabs):
        x = -1 - i
        script += f"DETECTOR({i + 0.5}, {0}, 0) rec[{x}] rec[{x - num_stabs}]\n"


    script += f"M {' '.join(map(str, data_qubits))}\n"

    for i in range(d):
        script += f'OBSERVABLE_INCLUDE({i}) rec[{-i -1}]\n'

    return stim.Circuit(script)



def estimate_msp(p,
                 q,
                 p2,
                 p_prep,
                 d,
                 rounds,
                 num_shots: int = 1_000):
    
    stim_circuit = get_ghz_circuit(d, p, q, p2, p_prep, rounds)

    model = stim_circuit.detector_error_model(decompose_errors=False, approximate_disjoint_errors=True)
    matching = pymatching.Matching.from_detector_error_model(model)
    sampler = stim_circuit.compile_detector_sampler()
    syndrome, actual_observables = sampler.sample(shots=num_shots, separate_observables=True)

    predicted_observables = matching.decode_batch(syndrome)

    errors = np.logical_xor(predicted_observables, actual_observables)

    num_err = abs(d - 2*np.sum(errors, axis = 1))

    rms_phase = np.sum(num_err**2)/num_shots

    return rms_phase


print('Generating data for average mean square phase....')


num_shots = 20000 # 20_000
from tqdm import tqdm
p = np.linspace(0.001,0.06,10)
p2 = p/5

d = np.array([7, 21 , 35, 41, 55], dtype = int)
L =  np.array(5*np.ceil(np.log(d)), dtype = int)
# L = np.copy(d) 

phases = []

for j in range(len(d)):
    phase = np.zeros_like(p,dtype = float)
    for i in tqdm(range(len(p))):
        phase[i] = estimate_msp(p[i], p[i], p2[i], p[i], d[j], L[j] ,num_shots)
    phases.append(phase)

np.savetxt('data/msp.txt', phases)

print('Data saved.')


print('Finding threshold...')


p = 0.06
dp = 0.001
p2 = p/5
d1 = 21; d2 = 42
num_shots = 20000


r1 = np.sqrt(estimate_msp(p, p, p2, p, d1, 5*np.ceil(np.log(d1)) ,num_shots)/d1**2)
r2 = np.sqrt(estimate_msp(p, p, p2, p, d2, 5*np.ceil(np.log(d1)) ,num_shots)/d2**2)

print('average mean square phase value at system size 21 and 42, respectively')

print(r1, r2)

while (r2-r1)>= 0:

    p += dp
    r1 = np.sqrt(estimate_msp(p, p, p2, p, d1, 5*np.ceil(np.log(d1)) ,num_shots)/d1**2)
    r2 = np.sqrt(estimate_msp(p, p, p2, p, d2, 5*np.ceil(np.log(d1)) ,num_shots)/d2**2)

    print(r1, r2)

print('Cross-over found at ', p)


print('Generating plot...')

fig, ax = plt.subplots(figsize = (5,3.5))

import seaborn as sns
# Use seaborn's crest palette
pal = sns.color_palette("rocket_r", n_colors=len(d))

for i in range(len(d)):

    phase = phases[i]
    
    ax.plot(
        p, np.sqrt(phase/d[i]**2),
        linestyle='-',
        linewidth=2,
        color=(*pal[i], 0.4),              # line color
        marker='o',
        markersize=6,
        markerfacecolor=(*pal[i],1),  # marker fill
        markeredgecolor=(*pal[i],1),
        # alpha=1,
        label=r'$n = $ ' + f'{d[i]}',
    )

ax.axvline(0.047, color='gray', linestyle='-.', linewidth=1.5, zorder = 0)


# Set the color of all spines to gray
for spine in ax.spines.values():
    spine.set_edgecolor('gray')

# plt.yticks([0,0.2,0.4,0.6],[0,0.2,0.4,0.6])

plt.tick_params(top=True, left=True, right=True, bottom=True,direction="in",axis='both', which='both', labelsize=15,color = 'k')
font2 = {'color':'black','size':18}
plt.xlabel("Physical error rate $p$",fontdict = font2)
plt.ylabel(r"$m_\mathrm{rms}/n$",fontdict = font2)
plt.legend(fontsize  = 15, ncol = 1, loc = 'best', frameon = False)
# plt.savefig('msp-correlated.svg',dpi = 500, bbox_inches='tight')
plt.savefig('data/msp-correlated.png',dpi = 500, bbox_inches='tight')
plt.show()