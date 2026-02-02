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
    
    av_msp = 0
    
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

        av_msp += pow(phase, 2)/num_shots
        
    return av_msp 



### State-preparation threshold 

print("\033[31mCaution\033[0m", ': The data is being generated for smaller number for shots for quick results. You can change its value to higher value on the script.')



print('Generating data for mean square phase...')

num_shots = 5_000

p = np.linspace(0.01,0.1,10)

d = np.array([7, 15, 21 , 35, 42], dtype = int)
L = np.array(5*np.ceil(np.log(d)), dtype = int)

phases = np.zeros((len(p), len(d)), dtype = float)


for i in tqdm(range(len(p))):
    for j in range(len(d)):
        phases[i,j] = estimate_mean_sq_phase(p[i],p[i],d[j],L[j],num_shots)


np.savetxt('data/mean-sq-phase-log-p_eq_q.txt', phases)

print('Done. Saved!')





### Uncomment for unequal error rates

print("Uncomment the following lines to generate data for unequal error rates. Skipping...")

# # p = np.linspace(0.01, 0.12, 12)
# p = np.linspace(0.01, 0.06, 10)
# d = np.array([7, 15, 21, 35, 42], dtype = int)
# L = np.array(5*np.ceil(np.log(d)), dtype = int)


# phases = np.zeros((len(L),len(p)), dtype = float)

# for j in range(len(L)):
#     for i in tqdm(range(len(p))):        
#         phases[j,i] = estimate_mean_sq_phase(p[i], 2*p[i], d[j], L[j], num_shots)


# print('code ran successfully.')

# np.savetxt('data/mean_sq_phases-measurement-p-twice.txt',phases)




# p = np.linspace(0.01, 0.12, 12)
# d = np.array([7, 15, 21, 35, 42], dtype = int)
# L = np.array(5*np.ceil(np.log(d)), dtype = int)


# phases = np.zeros((len(L),len(p)), dtype = float)

# for j in range(len(L)):
#     for i in tqdm(range(len(p))):        
#         phases[j,i] = estimate_mean_sq_phase(p[i], p[i]/2, d[j], L[j], num_shots)


# print('code ran successfully.')


# np.savetxt('data/mean_sq_phases-measurement-p-half.txt',phases)



print('Finding the threshold using finite system-size scaling...')


p = 0.06
dp = 0.001

d1 = 21; d2 = 42
num_shots = 5000


r1 = np.sqrt(estimate_mean_sq_phase(p,p,d1,int(5*np.ceil(np.log(d1))),num_shots)/d1**2)
r2 = np.sqrt(estimate_mean_sq_phase(p,p,d2,int(5*np.ceil(np.log(d2))),num_shots)/d2**2)

print('average mean square phase value at system size 21 and 42, respectively')

print(r1, r2)

while (r2-r1)>= 0:

    p += dp

    r1 = np.sqrt(estimate_mean_sq_phase(p,p,d1,int(5*np.ceil(np.log(d1))),num_shots)/d1**2)
    r2 = np.sqrt(estimate_mean_sq_phase(p,p,d2,int(5*np.ceil(np.log(d2))),num_shots)/d2**2)
    print(r1, r2)

print('Cross-over found at ', p)
  



print('Generating data for the mean square phase scaling below threshold...')


num_shots = 5000

p = [0.01, 0.02, 0.03, 0.04, 0.05]

d = np.array([10, 20, 50, 75, 100, 150],dtype = int)
L = np.array(5*np.ceil(np.log(d)), dtype = int)

phases = np.zeros((len(p), len(d)), dtype = float)

for i in tqdm(range(len(p))):
    for j in range(len(d)):
        phases[i,j] = estimate_mean_sq_phase(p[i],p[i],d[j],L[j],num_shots)


np.savetxt('data/mean_sq_phases_below_threshold.txt', phases)
print('Data saved.')



print('Generating plots...')



import matplotlib
matplotlib.use("Agg")  # or "Qt5Agg"

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

pal = ["003049","d62828","f77f00","fcbf49","eae2b7"]
pal = ['#' + i for i in pal]

marks = ['o', 's', '^', 'D', 'v', '>', '<', 'p', '*', 'h', 'H', '+', 'x']
import seaborn as sns

from matplotlib import font_manager
import matplotlib as mpl

num_shots = 20_000

p = np.linspace(0.01,0.1,10)

d = np.array([7, 15, 21 , 35, 42], dtype = int)

L = np.array(5*np.ceil(np.log(d)), dtype = int)


phases = np.loadtxt('data/mean-sq-phase-log-p_eq_q.txt')

fig, ax = plt.subplots(figsize = (5,4))

# Use seaborn's crest palette
pal = sns.color_palette("rocket_r", n_colors=len(d))

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Google Sans']


for i in range(len(d)):

    phase = phases[:,i]

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


ax.axvline(0.067, color='gray', linestyle='-.', linewidth=1.5, zorder = 0)

# Set the color of all spines to gray
for spine in ax.spines.values():
    spine.set_edgecolor('gray')

plt.tick_params(top=True, left=True, right=True, bottom=True,direction="in",axis='both', which='both', labelsize=15,color = 'k')
font2 = {'color':'black','size':18}
plt.xlabel("Error rate $\mathsf{p}$",fontdict = font2)
plt.ylabel(r"$m_\mathrm{rms}/n$",fontdict = font2)
plt.legend(fontsize  = 15,ncol = 1,loc = 'best',frameon = False)
# plt.savefig('mean-sq-phase.svg',dpi = 500, bbox_inches='tight')
plt.show()







num_shots = 5000

p = [0.01, 0.02, 0.03, 0.04, 0.05]

d = np.array([10, 20, 50, 75, 100, 150],dtype = int)
L = np.array(5*np.ceil(np.log(d)), dtype = int)

phases = np.loadtxt('data/mean_sq_phases_below_threshold.txt')

from scipy.optimize import curve_fit

def overx(x, a, b, c):
    return 1 - a - b/x**c

def linear(x, m):
    return m*x 


Ocean_sunset = ["#355070","#6d597a","#b56576","#e56b6f","#eaac8b"]


fig, ax = plt.subplots(figsize = (5,4))

n_col = 7
pal = sns.cubehelix_palette(as_cmap=False, n_colors = n_col)
blue = '#14213d'
# Plot with error bars

saturation = []

parameters = []

for i in range(len(p)):

    phase = phases[i,:]

    y = np.sqrt(phase/d**2)
    col = Ocean_sunset[i ] #purple[n_col - i - 2]
    ax.errorbar(d, y,
                color=col, 
                # label=r'Sampled logical error rate at $p = $ ' + f'{p}',
                marker=marks[0 % len(marks)],
                markersize=8,
                linestyle='',
                alpha=1,
                label = rf'$p =$ {p[i]}')

    gray = '#8d99ae'
    params, covariance = curve_fit(overx, d, y)

    x = np.arange(d[0],d[-1],0.01)
    ax.plot(x, overx(x, params[0], params[1], params[2]), color = col, linewidth = 3, alpha = 0.75)
    saturation.append(params[0])


ax.set_yscale("log")
ax.set_xscale("log")

# Set the color of all spines to gray
for spine in ax.spines.values():
    spine.set_edgecolor('gray')


plt.tick_params(top=True, left=True, right=True, bottom=True,direction="in",axis='both', which='both', labelsize=15,color = 'k')

font2 = { 'color':'black', 'size':18}

ax.set_xticks([10, 50, 100, 150],[10, 50, 100, 150])


ax.set_yscale("log")
ax.set_yticks([0.94, 0.96, 0.98, 1])
ax.set_yticklabels([0.94, 0.96, 0.98, 1])
ax.yaxis.set_minor_locator(plt.NullLocator())   # remove log minor ticks


plt.xlabel(rf"System size $n$", color = 'black', fontdict=font2)
plt.ylabel(r"$m_\mathrm{rms}/n$",fontdict = font2)
plt.show()



fig, ax = plt.subplots(figsize = (5,4))

blue = '#14213d'
# Plot with error bars

def square(x, a):
    return a*x**2



gray = '#adb5bd'
yellow = '#fcbf49'
params, covariance = curve_fit(square, p, saturation)
xi = np.arange(p[0], 0.05+0.0001, 0.0001)
ax.plot(xi, square(xi, params[0]), color = gray, linewidth = 4, zorder = 0, linestyle = '-.', alpha = 1)

for i in range(len(p)): 
    col = Ocean_sunset[i]
    ax.scatter(p[i], saturation[i], c = [col], s = 60, zorder = 1)


# Set the color of all spines to gray
for spine in ax.spines.values():
    spine.set_edgecolor('gray')


plt.tick_params(top=True, left=True, right=True, bottom=True,direction="in",axis='both', which='both', labelsize=15,color = 'k')
font1 = {'family':'serif', 'color':'black', 'size':22}
font2 = {'family':'serif', 'color':'black', 'size':18}
plt.xlabel(r"Error rate $\mathsf{p}$", color = 'black', fontdict=font2)
plt.ylabel(r"$\gamma$", fontdict=font1)
plt.show()



# Ucomment to plot unequal probability case - 


# num_shots = 5000


# p = np.linspace(0.01, 0.12, 12)

# d = np.array([7, 15, 21, 35, 42], dtype = int)
# L = np.array(5*np.ceil(np.log(d)), dtype = int)



# phases = np.loadtxt('data/mean_sq_phases-measurement-p-half.txt')

# fig, ax = plt.subplots(figsize = (3.5,5))

# # Use seaborn's crest palette
# pal = sns.color_palette("rocket_r", n_colors=len(d))

# mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['font.sans-serif'] = ['Google Sans']


# for i in range(len(d)):

#     phase = phases[i,:]

#     ax.plot(
#         p, np.sqrt(phase/d[i]**2),
#         linestyle='-',
#         linewidth=2,
#         color=(*pal[i], 0.4),              # line color
#         marker='o',
#         markersize=6,
#         markerfacecolor=(*pal[i],1),  # marker fill
#         markeredgecolor=(*pal[i],1),
#         # alpha=1,
#         label=r'$n = $ ' + f'{d[i]}',
#     )

# ax.axvline(0.097, color='gray', linestyle='-.', linewidth=1.5, zorder = 0)



# # Set the color of all spines to gray
# for spine in ax.spines.values():
#     spine.set_edgecolor('gray')

# plt.xticks([0.01,0.05, 0.1])

# plt.yticks([0.8, 0.9, 1])


# plt.tick_params(top=True, left=True, right=True, bottom=True,direction="in",axis='both', which='both', labelsize=15,color = 'k')
# font2 = {'color':'black','size':18}
# plt.xlabel(r"Error rate $\mathsf{p}$",fontdict = font2)
# plt.ylabel(r"$m_\mathrm{rms}/n$",fontdict = font2)
# plt.legend(fontsize  = 15,ncol = 1,loc = 'best',frameon = False)
# plt.show()





# num_shots = 5000


# # p = np.linspace(0.01, 0.12, 12)
# p = np.linspace(0.01, 0.06, 10)
# d = np.array([7, 15, 21, 35, 42], dtype = int)
# L = np.array(5*np.ceil(np.log(d)), dtype = int)



# phases = np.loadtxt('data/mean_sq_phases-measurement-p-twice.txt')

# fig, ax = plt.subplots(figsize = (3.5, 5))

# # Use seaborn's crest palette
# pal = sns.color_palette("rocket_r", n_colors=len(d))

# mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['font.sans-serif'] = ['Google Sans']


# for i in range(len(d)):

#     phase = phases[i,:]


#     ax.plot(
#         p, np.sqrt(phase/d[i]**2),
#         linestyle='-',
#         linewidth=2,
#         color=(*pal[i], 0.4),              # line color
#         marker='o',
#         markersize=6,
#         markerfacecolor=(*pal[i],1),  # marker fill
#         markeredgecolor=(*pal[i],1),
#         # alpha=1,
#         label=r'$n = $ ' + f'{d[i]}',
#     )

# ax.axvline(0.044, color='gray', linestyle='-.', linewidth=1.5, zorder = 0)


# # Set the color of all spines to gray
# for spine in ax.spines.values():
#     spine.set_edgecolor('gray')


# plt.xticks([0.01,0.03, 0.05])

# plt.yticks([0.8, 0.9, 1])


# plt.tick_params(top=True, left=True, right=True, bottom=True,direction="in",axis='both', which='both', labelsize=15,color = 'k')
# font2 = {'color':'black','size':18}
# plt.xlabel(r"Error rate $\mathsf{p}$",fontdict = font2)
# plt.ylabel(r"$m_\mathrm{rms}/n$",fontdict = font2)
# plt.legend(fontsize  = 15,ncol = 1,loc = 'best',frameon = False)
# plt.show()
