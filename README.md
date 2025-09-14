chickenpox_sindy.py implements the SINDy algorithm to model the observed rates of Chickenpox infection rates in Onatrio, Canada, using differential equations through sparse regression.
Result: Derives differential equations that govern the spread of Chickenpox using the SIR model.

Set up:
Datasets considered to create the SINDy models for Chickenpox and Rubella:
Weekly data of Chickenpox infection incidence in Ontario from 1946 to 1967
Weekly data of Rubella infection incidence in Ontario from 1946 to 1960
Population data of Ontario in 1931, 1941, 1951, 1956, 1961, 1966 and 1971, which was subjected to linear interpolation to generate weekly data from 1946 to 1967
Monthly data of births in Ontario from 1946 to 1967, which was subjected to linear interpolation to generate weekly data for the same period.
Application of a Savitzky-Golay filter of order 3 with a window length of 19 to reduce the possibility of SINDy overfitting data 


Recreating susceptibility time series:

Using the global regression technique, a time series of susceptibility data was recreated based on the following components:
St = S0 + Zt
S0: initial size of population assumed to be susceptible (~0.1 for Chickenpox and Rubella)St = proportion of population that is susceptible at time t
Zt = Y(end) -  alpha*X(end)
Yt = ∑ Bi,i+1                                 Bi,i+1 represents the number of new births in the interval from i to i+1
Xt = ∑ Ci,i+1                                Ci,i+1 represents the number of new cases in the interval from i to i+1
alpha = reporting rate
The Susceptibility time series obtained by these components was then normalized with respect to corresponding population data

<img width="1218" height="756" alt="image" src="https://github.com/user-attachments/assets/c836f61d-7b93-4231-a6d1-d5a8ebf6c573" />
<img width="745" height="498" alt="image" src="https://github.com/user-attachments/assets/a0bf4e1f-548d-4e6a-9cc8-3d2fcc782078" />



Recreating transmission rate time series:

A time series for the transmission rate was created using the FC method. (Equation 7 of https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008124)
βt ≈ (Zt, t+Δt) / (St It Δt)βt represents the transmission rate at time t
Z t, t+Δt   represents the number of infections that occur during the interval (t, t+Δt)
St represents the number of susceptible people at time t
It represents the number of infected people at time t

<img width="1184" height="765" alt="image" src="https://github.com/user-attachments/assets/0e94e743-7572-4b7a-9acf-bb299adc63c1" />

<img width="887" height="665" alt="image" src="https://github.com/user-attachments/assets/b364a2ae-2061-4863-98ba-d5207ceb162f" />


(S)' = -0.426 + -3.013 S + 3664.936 I + 0.132 β + 138.714 S^2 + -14078.294 S I + -1.615 S β + 494436.607 I^2 + -482.233 I β + -0.005 β^2 + -681.426 S^3 + -55393.210 S^2 I + 4.212 S^2 β + 1620.533 S I β + 0.037 S β^2 + -24834.715 I^2 β + 12.238 I β^2
(I)' = 0.039 + -0.713 S + -103.092 I + -0.005 β + 4.140 S^2 + 390.311 S I + 0.061 S β + 8.631 I β + -8.082 S^3 + -0.174 S^2 β + 1.013 S I β + -0.001 S β^2 + -0.235 I β^2
<img width="3991" height="170" alt="image" src="https://github.com/user-attachments/assets/c4e4ced3-f209-4365-82b1-a50dfe3b3890" />






