# secure-asynchronous-state-estimation
MATLAB simulation code for secure state estimation in LTI systems with asynchronous and non-periodic measurements

### Requirements
- MATLAB R2023 or later
- YALMIP toolbox (https://yalmip.github.io/tutorial/installation/)
- mosek solver (https://www.bing.com/search?q=mosek&qs=n&form=QBRE&sp=-1&lq=0&pq=mosek&sc=13-5&sk=&cvid=91562ECDC09B48479E2A648C52A8ED5C)

### How to Run
- Open the simulation file: TAC_Secure_Async_State_Estimation.m
- edit ns, the number of time steps, you want to simulate
- edit gammaLs, regularization parameter of the secure least square optimization, to trade-off the estimation performance 
