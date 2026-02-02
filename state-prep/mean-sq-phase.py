from essentials import *

def estimate_mean_sq_phase(qubit_error: float, measurement_error: float, distance: int, num_rounds: int = None, num_shots: int = 1_000):

    if num_rounds is None:
        num_rounds = 2*distance

    code = RepetitionCode(distance, num_rounds - 1, qubit_error, measurement_error)
    nodes = code.nodes

    initializer = code.get_initializer()
    solver = fb.SolverSerial(initializer)


    H = repetition_stabilisers(distance)

    matching = Matching(H)
    
    average_phase = 0
    
    for _ in range(num_shots):

        qubit_hist = np.tile(np.random.randint(2,size = (distance,1)),(1,num_rounds)) # create noise-less qubit history
        qubit_noise = (np.random.rand(distance, num_rounds) < qubit_error).astype(np.uint8) 
        qubit_noise[:,0] = 0 # first column noise doesn't mean anything 
        noise_cumulative = (np.cumsum(qubit_noise, 1) % 2).astype(np.uint8)
        noisy_qubit_hist = (qubit_hist + noise_cumulative)%2

        
        S = H@noisy_qubit_hist % 2

        syndrome_noise = (np.random.rand(distance-1, num_rounds) < measurement_error).astype(np.uint8)

        noisy_syndrome = (S + syndrome_noise) % 2

        syndrome_change = noisy_syndrome[:,1:] ^ noisy_syndrome[:,:-1] # (noisy_syndrome[:,1:] - noisy_syndrome[:,0:-1])%2 # add detectors

        Se = np.pad(syndrome_change, pad_width=1, mode='constant', constant_values=0)


        # get the indices where the value is 1

        row_indices, col_indices = np.where(Se == 1)

        defects = nodes[row_indices, col_indices].tolist()

        # defects = []

        # for i in range(len(row_indices)):
        #     defects.append(nodes[row_indices[i],col_indices[i]])
            
        syndrome = fb.SyndromePattern(defect_vertices=defects)

        solver.solve(syndrome)
        subgraph = solver.subgraph()

        spacelike_labels, timelike_labels = code._edge_labels()


        predicted_syndrome_change = np.zeros(distance - 1,dtype = int) 

        for i in range(distance - 1):
            if timelike_labels[i,-1] in subgraph:
                predicted_syndrome_change[i] = 1

        predicted_syndrome = (noisy_syndrome[:,-1] + predicted_syndrome_change)%2
        
    # -------------------------

        actual_syndrome = S[:,-1]
        
        
        solver.clear()

        correct_correction = np.array(matching.decode(actual_syndrome), dtype = float)
        predicted_correction = np.array(matching.decode(predicted_syndrome), dtype = float)

        error = np.abs(correct_correction - predicted_correction)
        
        phase = np.abs(np.count_nonzero(error == 0) - np.count_nonzero(error == 1))

        average_phase += pow(phase, 2)/num_shots
        
    return average_phase 



num_shots = 20_000

# p = np.linspace(0.01,0.1,10)

# d = np.array([7, 15, 21 , 35, 42], dtype = int)
# L = np.array(5*np.ceil(np.log(d)), dtype = int)

# phases = np.zeros((len(p), len(d)), dtype = float)


# for i in tqdm(range(len(p))):
#     for j in range(len(d)):
#         phases[i,j] = estimate_mean_sq_phase(p[i],p[i],d[j],L[j],num_shots)


# np.savetxt('mean-sq-phase-log-p_eq_q.txt', phases)


# # p = np.linspace(0.01, 0.12, 12)
# p = np.linspace(0.01, 0.06, 10)
# d = np.array([7, 15, 21, 35, 42], dtype = int)
# L = np.array(5*np.ceil(np.log(d)), dtype = int)


# phases = np.zeros((len(L),len(p)), dtype = float)

# for j in range(len(L)):
#     for i in tqdm(range(len(p))):        
#         phases[j,i] = estimate_mean_sq_phase(p[i], 2*p[i], d[j], L[j], num_shots)


# print('code ran successfully.')

# # print(phases)

# np.savetxt('mean_sq_phases-measurement-p-twice.txt',phases)





p = 0.03
dp = 0.001

d1 = 21; d2 = 42
num_shots = 20_000


r1 = np.sqrt(estimate_mean_sq_phase(p,2*p,d1,int(5*np.ceil(np.log(d1))),num_shots)/d1**2)
r2 = np.sqrt(estimate_mean_sq_phase(p,2*p,d2,int(5*np.ceil(np.log(d2))),num_shots)/d2**2)

while (r2-r1)>= 0:

    p += dp

    r1 = np.sqrt(estimate_mean_sq_phase(p,2*p,d1,int(5*np.ceil(np.log(d1))),num_shots)/d1**2)
    r2 = np.sqrt(estimate_mean_sq_phase(p,2*p,d2,int(5*np.ceil(np.log(d2))),num_shots)/d2**2)

    print(p)
  

# num_shots = 20_000

# p = [0.01, 0.02, 0.03, 0.04, 0.05]

# d = np.array([10, 20, 50, 75, 100, 150],dtype = int)
# L = np.array(5*np.ceil(np.log(d)), dtype = int)

# phases = np.zeros((len(p), len(d)), dtype = float)

# for i in tqdm(range(len(p))):
#     for j in range(len(d)):
#         phases[i,j] = estimate_mean_sq_phase(p[i],p[i],d[j],L[j],num_shots)


# np.savetxt('mean_sq_phases_below_threshold.txt', phases)


