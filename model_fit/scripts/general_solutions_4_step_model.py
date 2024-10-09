# analytical solutions to the general, rescaled 4-step model

import numpy as np

# --- DEFINE FUNCTIONS FOR THE SOULTIONS TO THE 4-STEP MODEL ----------------
# for the 3-step model simply the first three funtions need to be used
# for single-exon genes the step1 is used for mature mRNA (so one less step is fitted)

# STEP 1:
def step1(t, r1):
    x = 1.-np.exp(-t*r1)
    return x


# STEP 2:
def step2(t, r1, r2):
    x = (((1.-(np.exp(-t*r2)))*r1)+((-1.+(np.exp(-t*r1)))*r2))/(r1-r2)
    return x


# STEP 3:
def step3(t, r1, r2, r3):
    x = ((1.+((((np.exp(-t*r2))-(np.exp(-t*r3)))*r3)/(r2-r3)))-
         ((((r2*(r3*(((np.exp(-t*r3))*(r1-r2))+
                     (((np.exp(-t*r1))*(r2-r3))+((np.exp(-t*r2))*(r3-r1))))))/
                     (r2-r3))/(r1-r3))/(r1-r2)))-(np.exp(-t*r3))
    return x


# STEP 4:
def step4(t, r1, r2, r3, r4):
    x = (((1.+((((np.exp(-t*r3))-(np.exp(-t*r4)))*r4)/(r3-r4)))-
          (r2*(r3*(r4*((((((np.exp(-t*r2))/(r2-r4))/(r2-r3))/(r1-r2))+
                        (((((np.exp(((-t*r3))))/(r3-r4))/(r3-r2))/(r1-r3))+
                         ((((np.exp(-t*r4))/(r4-r3))/(r4-r2))/(r1-r4))))-
                         ((((np.exp(-t*r1))/(r1-r4))/(r1-r3))/(r1-r2)))))))-
                         ((((r3*(r4*(((np.exp(-t*r4))*(r2-r3))+
                                     (((np.exp(-t*r2))*(r3-r4))+
                                      ((np.exp(-t*r3))*(r4-r2))))))/
                                      (r3-r4))/(r2-r4))/(r2-r3)))-(np.exp(-t*r4))
    return x


# --- function that accounts for time delays until T2C conversions are seen ----
def modified_time(t, elongation_rate, distance_to_end, overall_delay):
    t_mod = t - distance_to_end / elongation_rate - overall_delay
    t_mod[ t_mod < 0 ] = 0.
    return t_mod


# --- simpler modified time that accounts only for overall delay until labeling is seen ----
def delayed_time(t, overall_delay):
    t_mod = t - overall_delay
    t_mod[ t_mod < 0 ] = 0.
    return t_mod

