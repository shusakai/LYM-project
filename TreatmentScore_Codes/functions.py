### import modules
import numpy as np
import pandas as pd
from scipy import constants

### cell radius (reused from Shilts et al. Nature 2022)
cell_rad = {'B.memory':3, 'B.naive':3, 'B.plasma':3, 'Treg.memory':4, 'Treg.naive':4, 'T4.naive':4, 'T4.CM':4, 'T4.EM':4, 'T4.EMRA':4, 'Th1':4, 'Th17':4, 'Th2':4, 'T8.naive':4, 'T8.CM':4, 'T8.EM':4, 'T8.EMRA':4, 'mDC':5, 'pDC':5, 'Mono.classical':8, 'Mono.intermediate':8, 'Mono.nonclassical':8, 'NK.bright':8, 'NK.dim':8}

########################################################################
# These three functions are reused from Shilts et al. Nature 2022.
def calculate_SA(cell_radius):
    ''' input cell radius and outputs surface area in the same units squared'''
    return 4 * constants.pi * (cell_radius)**2

def calculate_3Dto2D(Kd_3D, height = 0.010):
    ''' Converts 3D Kd in Molar units to 2D Kd in units moelcules/um2 (height is in um units)  \n
        uses constant height factor to get rough estimates since very few are empirically measured'''
    return constants.Avogadro * 1E-15 * height * Kd_3D * 1E-9
    # The unit of 3D Kd in the reference is nM, hence 1E-9 is needed
    
def bound_concentration(conc_A_total, conc_B_total, Kd_2D):
    '''
    Calculates amount of bound A-B given the total amount of A and B and their 2D affinity \n
     dervied from the quadratic solution the system of equations \n
     Kd = Rec * Lig / Bound \n
     Rec = Rec_tot - Bound  \n
     Lig = Lig_tot - Bound  \n
    Where all concentrations and Kd are in the same units (e.g. molecules/um2)
    Beware this equation does not handle case if 0 is input as a concentration or Kd (because smaller numbers = stronger binding so zero is perfect conversion to bound state)
    '''
    conc_AB = ((Kd_2D + conc_A_total + conc_B_total) -
               ((Kd_2D + conc_A_total + conc_B_total)**2 -
               (4 * conc_A_total * conc_B_total))**(1/2)) / 2
    return conc_AB
########################################################################


# Return [total score, sum of scores for target interaction]
def interaction_score(c1, c2, data1, data2, df1, antag):
    l1 = list(data1.index)
    l2 = list(data2.index)
    l_gene = list( set(l1) & set(l2) )
    st1 = c1.split('_')[0]
    st2 = c2.split('_')[0]
    cell_SA1 = calculate_SA(cell_rad[st1])
    cell_SA2 = calculate_SA(cell_rad[st2])
    df_interaction = df1
    x1 = np.zeros(2)
    
    for row in df_interaction.itertuples(name=None):
        #print(row[1], row[2], row[3])   # molecule1, molecule2, KD
        if (row[1] in l_gene) and (row[2] in l_gene):
            if data1.at[row[1], st1] and data1.at[row[2], st2]:
                K_D = calculate_3Dto2D(row[3], height = 0.010)
                conc_1_total = data2.at[row[1], c1] / cell_SA1  # molecule1 on c1
                conc_2_total = data2.at[row[2], c2] / cell_SA2  # molecule2 on c2
                conc = bound_concentration(conc_1_total, conc_2_total, K_D)
                x1[0] = x1[0] + conc
                if (row[1] in antag) or (row[2] in antag):
                    x1[1] = x1[1] + conc
            
            if data1.at[row[1], st2] and data1.at[row[2], st1]:
                K_D = calculate_3Dto2D(row[3], height = 0.010)
                conc_1_total = data2.at[row[1], c2] / cell_SA2  # molecule1 on c2
                conc_2_total = data2.at[row[2], c1] / cell_SA1  # molecule2 on c1
                conc = bound_concentration(conc_1_total, conc_2_total, K_D)
                x1[0] = x1[0] + conc
                if (row[1] in antag) | (row[2] in antag):
                    x1[1] = x1[1] + conc
                
    if c1 == c2:
        # interaction between same cell type will return twiced score
        return x1/2.0
    else:
        return x1

