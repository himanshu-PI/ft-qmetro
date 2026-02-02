from essentials import *

from math import comb

from scipy.optimize import curve_fit

# Define the quadratic function to fit (a * x^2)
def overx(x, a, b, c):
    return 1 - a - b/x**c

def linear(x, m):
    return m*x 

def errorX(a, q):
    
    threshold = a // 2
    prob = sum(comb(a, k) * (q**k) * ((1-q)**(a-k)) for k in range(threshold + 1))

    return prob


def errorLogX(a, q, n):
    """
    Logical X decoding error for n independent physical X majority-vote blocks.

    a : number of repetitions per physical measurement
    q : physical measurement error probability
    n : number of independent X measurements

    Returns: probability all n are correct (success probability).
    """

    p_succ = errorX(a, q)     # success of a single majority-vote block
    return np.min([1 - p_succ**n , 0.5])         # all n must succeed


def FIqec(errMeas, errPrep, errQubit , msp, n):

    x1 = 1 - 2*errMeas
    x2 = 1 - 2*errPrep
    x3 = 1 - 2*errQubit
    p = errQubit

    y1 = (x1 * x2 * x3)*((x3**2)*(msp/n) + 4*p*(1-p))

    return y1**2


# Define the quadratic function to fit (a * x^2)
def overx(x, a, b, c):
    return 1 - a - b/x**c

def linear(x, m):
    return m*x 



print("\033[31mCaution\033[0m", ': You must run mean-sq-phase.py to generate require data before running this code.')


print('Calculation FI for our protocol below threshold...')

const = 5

p = 0.04

p_prep = 0.01
q = 0.01

qx = p_prep*(1-q) + q*(1-p_prep)


n = 2**(np.arange(3, 9, 1))

msp = np.loadtxt(f'data/msp-p-{p}.txt')


errX = np.zeros_like(n, dtype = float)

for i in range(len(n)):
    rounds = int(const*np.log(n[i]))
    errX[i] = errorLogX(rounds, qx, n[i])


phases = np.loadtxt('data/mean_sq_phases_below_threshold.txt')
i = 1
phase = phases[i,:]
d = np.array([10, 20, 50, 75, 100, 150],dtype = int)
y = np.sqrt(phase/d**2)
params, covariance = curve_fit(overx, d, y)

msp_fit = pow(overx(n, params[0], params[1], params[2]), 2)* pow(n, 2)


fi_exact = np.zeros_like(errX)
fi_fit = np.zeros_like(errX)


for i in range(len(n)):

    fi_exact[i] = FIqec(errX[i], errX[i], p_prep, msp[i], n[i])
    fi_fit[i]  = FIqec(errX[i], errX[i], p_prep, msp_fit[i], n[i])



fig, ax = plt.subplots(figsize = (6,4))

# Use seaborn's crest palette
pal = sns.color_palette("crest", n_colors=4)

sql = n
hl = n**2

pal = ["#003049","#d62828","#f77f00","#fcbf49","#eae2b7"]

orange = '#e07a5f'
black = '#000000'

green = '#81b29a'
gray = '#adb5bd'
ax.plot(n, hl , linestyle = '-', color = orange, label = rf'$n^2$', linewidth = 2.5)
ax.plot(n, (1-2*p)*sql, linestyle = '-', color = pal[0], label = r'$(1-2q)n$', linewidth = 2.5)
ax.plot(n, (n/np.exp(1))*(1/np.log(1/(1-2*p)**2)), linestyle = '-', color = green, label = r'$\mathcal{F}_\mathrm{cl}^{(q)}$', linewidth = 2.5)


yellow = '#f2cc8f'

ax.plot(
    n, fi_exact,
    linestyle='-.',
    linewidth=2.5,
    color=yellow,              # line color
    marker='o',
    markersize=6,
    markerfacecolor='black',  # marker fill
    markeredgecolor='black',
    alpha=1,
    label=r'$\mathcal{F}_\mathrm{cl}^{(\mathrm{QEC})}$'
)


ax.set_yscale("log")
ax.set_xscale('log')

# Set the color of all spines to gray
for spine in ax.spines.values():
    spine.set_edgecolor('gray')

plt.tick_params(top=True, left=True, right=True, bottom=True,direction="in",axis='both', which='both', labelsize=15,color = 'k')
font1 = {'family':'serif', 'color':'black', 'size':12}
font2 = {'family':'serif','color':'black','size':18}
plt.xlabel(r"$n$",fontdict = font2)
plt.ylabel(r"$\mathcal{F}_\mathrm{cl}$",fontdict = font2)
plt.legend(fontsize  = 14,ncol = 1,loc = 'best',frameon = False)

ax.set_yscale("log")
plt.show()









const = 5

p = 0.1

n = 2**(np.arange(3, 9, 1))

sp = np.loadtxt(f'data/msp-p-{p}.txt')

errX = np.zeros_like(n, dtype = float)

for i in range(len(n)):
    rounds = int(const*np.log(n[i]))
    errX[i] = errorLogX(rounds, p, n[i])




fi_exact = np.zeros_like(errX)

for i in range(len(n)):

    fi_exact[i] = FIqec(errX[i], errX[i], p, msp[i], n[i])


fig, ax = plt.subplots(figsize = (6,4))

# Use seaborn's crest palette
pal = sns.color_palette("crest", n_colors=4)

sql = n
hl = n**2

pal = ["#003049","#d62828","#f77f00","#fcbf49","#eae2b7"]

orange = '#e07a5f'
black = '#000000'

green = '#81b29a'
gray = '#adb5bd'
ax.plot(n, hl , linestyle = '-', color = orange, label = rf'$N^2$', linewidth = 2)
ax.plot(n, (1-2*p)*sql, linestyle = '-', color = pal[0], label = r'$(1-2q)N$', linewidth = 2)
ax.plot(n, (n/np.exp(1))*(1/np.log(1/(1-2*p)**2)), linestyle = '-', color = green, label = r'$\mathcal{F}_\mathrm{cl}^{(q)}$', linewidth = 2)

# ax.plot(na, noisy_fi, linestyle = '-', color = pal[3], label = rf'$\eta^2 n^2$')

# ax.plot(n, rep_qubit**2, linestyle = '-', color = gray, label = r'$a$')


yellow = '#f2cc8f'

ax.plot(n, fi_exact,
            color='black', 
            marker= 'o',
            markersize=6,
            linestyle='',
            alpha=1,
            label=r'$\mathcal{F}_\mathrm{cl}^{(\mathrm{QEC})}$')



ax.set_yscale("log")
ax.set_xscale('log')

# Set the color of all spines to gray
for spine in ax.spines.values():
    spine.set_edgecolor('gray')

plt.tick_params(top=True, left=True, right=True, bottom=True,direction="in",axis='both', which='both', labelsize=15,color = 'k')
font1 = {'family':'serif', 'color':'black', 'size':12}
font2 = {'family':'serif','color':'black','size':18}
plt.xlabel(r"$N$",fontdict = font2)
plt.ylabel(r"$\mathcal{F}_\mathrm{cl}$",fontdict = font2)
plt.legend(fontsize  = 14,ncol = 1,loc = 'best',frameon = False)
ax.set_yscale("log")
plt.show()