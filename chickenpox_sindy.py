import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pysindy as ps
from pysindy import SINDy
from pysindy.feature_library import PolynomialLibrary
from pysindy.optimizers import STLSQ
import pandas as pd
from scipy.integrate import odeint
from scipy.stats import norm
import itertools
from scipy.signal import savgol_filter
from sklearn.metrics import mean_squared_error
import random


# Load your datasets
birth_data = np.loadtxt(r"Ontario_Birth_Data_M.txt")
chickenpox_data = np.loadtxt(r"OntarioChickenWeekly46_67.txt")
rubella_data = np.loadtxt(r"OntarioRubellaWeekly46_67.txt")
#pop_data = np.loadtxt(r"C:\Users\mauli\Documents\ILMEE\data for infectious disease\pop_ontario_47_67.txt")
pop_data = np.loadtxt(r"Ontario_Demographics_Measles.txt")
susc_chick = np.loadtxt(r"susc_chick.txt")

pop_time = pop_data[:,0]
pop = pop_data[:, 1]

susc_time = susc_chick[:,0]
susc_cases = susc_chick[:,1]

birth_time = birth_data[:, 0]
births = birth_data[:, 1]/12.8

chickenpox_time = chickenpox_data[:, 0][0:730]

pop_interp_func = interp1d(pop_time, pop, kind='linear')
interp_pop = pop_interp_func(chickenpox_time)
abs_chicken_cases = chickenpox_data[:, 1]
chickenpox_cases = chickenpox_data[:, 1][0:730]/interp_pop


rubella_time = rubella_data[:,0]
rubella_cases = rubella_data[:,1]


rub_interp_func = interp1d(rubella_time, rubella_cases, kind='linear')
rub_interp = rub_interp_func(chickenpox_time)
chickenpox_cases = rub_interp/interp_pop


chickenpox_cases = savgol_filter(chickenpox_cases, window_length = 19, polyorder=3)





transmission = np.zeros_like(chickenpox_cases)

#performing interpolation on births (which has data for every 4 months) to get data for every week corresponding to chickenpox_time)
births_interp_func = interp1d(birth_time, births, kind='linear')
interp_births = births_interp_func(chickenpox_time)

pop_interp_func = interp1d(pop_time, pop, kind='linear')
interp_pop = pop_interp_func(chickenpox_time)

susc_interp_func = interp1d(susc_time, susc_cases, kind='linear')
interp_susc = susc_interp_func(chickenpox_time)
#print(interp_susc[0])



years = np.array([
    1976, 1971, 1966, 1961, 1956, 1951, 1941, 1931, 1921, 1911, 1901,
    1891, 1881, 1871, 1861, 1851
])

ppl = np.array([
    8264465, 7703106, 6960870, 6236092, 5404933, 4597542, 3787655, 3431683,
    2933662, 2527292, 2182947, 2114321, 1926922, 1620851, 1396091, 952004
])

# Sorting the data as years are in descending order
sorted_indices = np.argsort(years)
years = years[sorted_indices]
ppl = ppl[sorted_indices]


interp_func = interp1d(years, ppl, kind='linear')

weekly_populations = interp_func(chickenpox_time)



def SuscRec_FGlocal(C, B, fac):
    # Calculate cumulative sums
    Y = np.cumsum(B)
    X = np.cumsum(C)
    
    P = np.polyfit(X, Y, 1)
    
    # Extract the slope coefficient (representing alpha)
    alpha = P[0]
    
    # Calculate the weighted difference
    Zt = Y[-1] - alpha * X[-1]

    return Zt, alpha
def SuscRec_FG_wrapper(C, B, fac, init):
    N = len(C)
    S0 = init

    # Initialize arrays
    S = np.zeros(N)
    alphaV = np.zeros(N)
    Z = np.zeros(N)

    # Set initial values
    S[0] = S0
    alphaV[0] = alpha = 8.5  # Example initial value for alpha
    Z[0] = 0
    # Iterate to compute Z and alphaV
    for i in range(1, N):
        Z[i], alphaV[i] = SuscRec_FGlocal(C[:i+1], B[:i+1], fac)
        S[i] = S0 + Z[i]

    # Normalize susceptible population
    S_normalized = S/interp_pop  # Assuming P is available
    
    return S_normalized




def rescale(array):
    scaled_cases = (array - np.min(array)) / (np.max(array) - np.min(array))
    rescaled_cases = (scaled_cases * (max(chickenpox_cases) - min(chickenpox_cases))) + min(chickenpox_cases)
    return rescaled_cases

def rescale_s(array):
    scaled_cases = (array - np.min(array)) / (np.max(array) - np.min(array))
    rescaled_cases = (scaled_cases * (max(S_normalized) - min(S_normalized))) + min(S_normalized)
    return rescaled_cases



def compute_aic(y_true, y_pred, num_params):
    resid = y_true[len(y_pred)//2:-1] - rescale(y_pred)[len(y_pred)//2:-1]
    sse = np.sum(resid ** 2)#/len(y_true)
    aic = 2 * num_params + len(y_true)//2 * np.log(2*sse / len(y_true))
    return sse
    #return aic


S0_values= range(0,round(interp_pop[0]),10000)

lambda_values = np.linspace(0.0000001, 0.1, 100000)
lambda_values = [0.0000001]

best_aic = float('inf')
best_model = None
best_params = (None, None)
aic_list ={}

window_s = random.randint(10,30)
window_t = random.randint(10,30)
poly_s = random.randint(1, 9)
poly_t = random.randint(1, 9)

for lam in [0.001]:#lambda_values:
    for init in [740000]:#]range(720000, 750000,1000):

        S_normalized = SuscRec_FG_wrapper(chickenpox_cases, interp_births, 0.45, init)

        for i in range(1, len(chickenpox_cases)):
            transmission[i] = 1*chickenpox_cases[i]/(S_normalized[i]*chickenpox_cases[i-1])
        transmission[0]=transmission[1]

        S_normalized = savgol_filter(S_normalized, window_length = 19, polyorder=3)

        transmission = savgol_filter(transmission, window_length = 29, polyorder=3)


        df = pd.DataFrame({
            'Time': chickenpox_time,
            'I': chickenpox_cases,
            'S': S_normalized,
            'β': transmission
        })
        #print(df)

        X = df[['S', 'I', 'β']].values
        t = df['Time'].values

        differentiation_method = ps.FiniteDifference()
        feature_library = ps.PolynomialLibrary(degree=3)
        optimizer = ps.optimizers.STLSQ(threshold=lam)   #play around with lamda

        model = ps.SINDy(
            differentiation_method=differentiation_method,
            feature_library=feature_library,
            optimizer=optimizer
        )

        model.feature_names=['S', 'I', 'β']
        model.fit(X, t=t)
        model.print()

        prediction = model.predict(X)
        #print("YO: ",prediction[:,1])
        print(init, lam)
        aic = compute_aic(X[:, 1], prediction[:, 1], len(model.coefficients()))
        aic_list[str(aic)] = (model, init, lam)
        # if min(prediction[:,1])>=0:
        #     break


        plt.figure(figsize=(10, 6))
        plt.plot(t, chickenpox_cases, label='true', color='blue')
        m = np.mean(savgol_filter(prediction[:,0], window_length = 80, polyorder=3))
        plt.plot(t, prediction[:,1], label='Infection Proportion', color ="red", linestyle='--')
        plt.xlabel("Time (Years)")
        plt.ylabel("Suscpetibility Proportion")
        plt.grid(True)
        plt.show()

        #break
        if aic<best_aic:
            best_aic = aic
            best_model = model
            best_params = (init, lam)
        #break
# print(best_aic)
# best_model.print()
# print(best_params)


nums=[float(x) for x in aic_list.keys()]

print(window_s)
print(window_t)
print(poly_s)
print(poly_t)
print(S0_values)





def diff_eq(y, t, params):
    S, I, beta = y

    dS_dt = -0.439 + 9.928*S + 3859.602*I + 0.042*beta - 52.324*S**2 - 26304.214*S*I - 0.877*S*beta - 378.007*I*beta - 0.001*beta**2 + 4.161*S**2*beta + 1512.337*S*I*beta + 0.009*S*beta**2 + 7.816*I*beta**2
    dI_dt = 0.001 - 0.025*S - 0.001*S*beta + 0.027*S**2*beta + 0.006*I*beta**2
    db_dt = -1393.086 + 64654.618*S - 95855875.217*I + 22.509*beta - 766610.500*S**2 + 1250849402.556*S*I - 1886.089*S*beta + 3314626818.538*I**2 + 6934621.586*I*beta + 1.390*beta**2 + 2665371.244*S**3 - 4163484613.748*S**2*I + 15126.837*S**2*beta - 45004882.392*S*I*beta - 1.186*S*beta**2 + 109118424.498*I**2*beta - 128802.226*I*beta**2 - 0.005*beta**3
    # dS_dt = -0.031 + 0.827*S + 1122.879*I - 0.001*beta - 3.094*S**2 - 8550.959*S*I - 0.001*S*beta - 54.277*I*beta
    # dI_dt = -0.001 + 0.002*S + 0.001*S*beta + 0.043*I*beta
    # db_dt = 223.199 - 4037.583*S + 233550.522*I - 6.934*beta + 22333.875*S**2 - 8710480.966*S*I + 53.571*S*beta + 5466243684.555*I**2 - 37805.709*I*beta + 0.110*beta**2
    # dS_dt = -0.031 + 0.827*S + 1122.879*I - 0.001*beta - 3.094*S**2 - 8550.959*S*I - 0.001*S*beta - 54.277*I*beta
    # dI_dt = 0.004 - 0.107*S + 0.526*S**2 + 0.003*S*beta + 0.048*I*beta
    # db_dt = 223.199 - 4037.583*S + 233550.522*I - 6.934*beta + 22333.875*S**2 - 8710480.966*S*I + 53.571*S*beta + 5466243684.555*I**2 - 37805.709*I*beta + 0.110*beta**2

    # dS_dt = -0.357 + 8.449*S - 263.759*I + 0.038*beta - 45.414*S**2 + 40799.962*S*I - 0.847*S*beta - 167.660*I*beta - 0.001*beta**2 - 277349.312*S**2*I + 4.171*S**2*beta - 153.230*S*I*beta + 0.009*S*beta**2 + 5.074*I*beta**2
    # dI_dt = 0.001 - 0.547*S**2 + 0.005*S*beta - 1.252*I*beta + 3.967*S**3 - 0.016*S**2*beta + 0.108*I*beta**2
    # db_dt = -1767.779 + 71173.867*S - 93959037.117*I + 57.081*beta - 781780.270*S**2 + 1149902028.226*S*I - 2338.552*S*beta + 53297137006.967*I**2 + 6873629.411*I*beta + 0.437*beta**2 + 2541217.951*S**3 - 3400295577.727*S**2*I + 16173.715*S**2*beta - 450469067177.364*S*I**2 - 43274024.952*S*I*beta + 6.338*S*beta**2 - 787259988.562*I**2*beta - 128751.765*I*beta**2 + 0.001*beta**3
    return [dS_dt, dI_dt, db_dt]
