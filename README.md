Liu, Y. & Li, G. Sure joint screening for high dimensional Cox's proportional hazards model under the case-cohort design.

Description:

This program implements the SMPLE joint screening method of the paper.

Main function: SMPLE(X,Y,delta,xi,beta_ini,k)

X: Design matrix

Y: Time to event or censoring

delta: 0=censoring; 1 otherwise

xi: Indicator being selected into the subcohort, 1 selected, 0 otherwise

beta_ini: Initial value for IWHT algorithm

k: Screened model size

S: Active index set

beta_S: Parameter estimates of active variables

Output: Active index set and parameter estimates of active variables.
