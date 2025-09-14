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


My recreation:
<img width="1218" height="756" alt="image" src="https://github.com/user-attachments/assets/c836f61d-7b93-4231-a6d1-d5a8ebf6c573" />

Recreation by Jonathan Horrocks & Chris T. Bauch:

<img width="745" height="498" alt="image" src="https://github.com/user-attachments/assets/a0bf4e1f-548d-4e6a-9cc8-3d2fcc782078" />



Recreating transmission rate time series:

A time series for the transmission rate was created using the FC method. (Equation 7 of https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008124)
βt ≈ (Zt, t+Δt) / (St It Δt)βt represents the transmission rate at time t
Z t, t+Δt   represents the number of infections that occur during the interval (t, t+Δt)
St represents the number of susceptible people at time t
It represents the number of infected people at time t

My recreation:
<img width="1184" height="765" alt="image" src="https://github.com/user-attachments/assets/0e94e743-7572-4b7a-9acf-bb299adc63c1" />

Recreation by Jonathan Horrocks & Chris T. Bauch:
<img width="887" height="665" alt="image" src="https://github.com/user-attachments/assets/b364a2ae-2061-4863-98ba-d5207ceb162f" />

Discovered SINDy model for Chickenpox using I, S and β as features upon fitting a 3rd degree polynomial:

(S)' = -0.426 + -3.013 S + 3664.936 I + 0.132 β + 138.714 S^2 + -14078.294 S I + -1.615 S β + 494436.607 I^2 + -482.233 I β + -0.005 β^2 + -681.426 S^3 + -55393.210 S^2 I + 4.212 S^2 β + 1620.533 S I β + 0.037 S β^2 + -24834.715 I^2 β + 12.238 I β^2
(I)' = 0.039 + -0.713 S + -103.092 I + -0.005 β + 4.140 S^2 + 390.311 S I + 0.061 S β + 8.631 I β + -8.082 S^3 + -0.174 S^2 β + 1.013 S I β + -0.001 S β^2 + -0.235 I β^2
<img width="3991" height="170" alt="image" src="https://github.com/user-attachments/assets/c4e4ced3-f209-4365-82b1-a50dfe3b3890" />


Chickenpox model inference:

(I)’ has a high dependence on I and the interaction between S and I i.e. SI

(S)’ has a high dependence on terms involving I, especially interactive terms like SI, I^2, S^2 I , S I β , and I^2 β  (Highest on I^2 and S^2 I)


Comparison of Real Chickenpox Data and Rescaled SINDy data:

<img width="1723" height="281" alt="image" src="https://github.com/user-attachments/assets/2f151faf-f768-4e28-9316-f3ca4ea2ea60" />

The trends depicted by the SINDy model are consistent with the peaks, ordering, and the cyclic nature of real infection data, which suggests that SINDy understand the dynamics of Chickenpox infections between 1946 and 1967 in Ontario


Comparison of Real Susceptibility Data and Rescaled SINDy data:
<img width="1666" height="358" alt="image" src="https://github.com/user-attachments/assets/678af324-1a5c-4036-8ac7-79dd54066f88" />


Recreated Susceptibility time series of Rubella:

<img width="1090" height="662" alt="image" src="https://github.com/user-attachments/assets/57e3a645-5269-4003-81dd-6e896ad223e7" />

Comparison with Jonathan Horrocks & Chris T. Bauch:


<img width="887" height="637" alt="image" src="https://github.com/user-attachments/assets/49b8dc12-c045-48b2-8e58-68fc1e2f1b6d" />



Comparison of recreated trend of average transmission rate with that of Jonathan Horrocks & Chris T. Bauch:
<img width="888" height="554" alt="image" src="https://github.com/user-attachments/assets/a5c443df-5a4f-4302-9e1d-942420284b31" />

<img width="840" height="672" alt="image" src="https://github.com/user-attachments/assets/738c5e34-6b83-4df6-9477-69489a5726af" />



Discovered SINDy model for Rubella using I, S and β as features upon fitting a 3rd degree polynomial:

(S)' = 17.399 + -290.357 S + -1320.308 I + 0.251 β + 1560.067 S^2 + 1.845 S β + -0.094 β^2 + -2673.927 S^3 + -17.882 S^2 β + 197.304 S I β + 0.475 S β^2 + 11.279 I β^2 + 0.001 β^3
(I)' = 0.014 I β^2
<img width="3576" height="164" alt="image" src="https://github.com/user-attachments/assets/5459444a-8ea6-4b86-be48-861178e4723f" />



Rubella model inferences:

(I)’ has a dependence on I β^2
(S)’ has the highest dependence on S, S^2, S^3, and S I β.



Comparison of Real Rubella Data and Rescaled SINDy data

<img width="1136" height="662" alt="image" src="https://github.com/user-attachments/assets/104db247-0966-475b-906c-46c5a66cf771" />


The trends depicted by the SINDy model are consistent with the peaks, ordering, and the cyclic nature of real infection data, which suggests that SINDy understand the dynamics of Rubella infections between 1946 and 1967 in Ontario


Comparison of Real Susceptibility Data and Rescaled SINDy data using a Savitzky-Golay filter:


<img width="1120" height="662" alt="image" src="https://github.com/user-attachments/assets/2cc1de4c-6ecb-4384-9f58-79b0c65caee7" />


Comparison of result with that of Jonathan Horrocks & Chris T. Bauch:


<img width="858" height="500" alt="image" src="https://github.com/user-attachments/assets/1e6a9ec6-ed8a-40b9-8dda-95b24ca54360" />

<img width="1128" height="467" alt="image" src="https://github.com/user-attachments/assets/f79d5463-86ea-40bc-a7fa-50e81cc8552c" />




