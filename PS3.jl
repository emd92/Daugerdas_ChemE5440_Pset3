cd("C:\\Users\\ellie\\OneDrive\\Desktop\\PS3")

#Problem Set 3
#Part A: Stoichiometric Matrix, S

S= [-1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; -1 0 0 1 2 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 -1 0 -2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0; 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 -1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; -1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0; 0 0 -1 0 4 -4 0 0 0 0 0 0 0 1 0 0 0 0 0 0 -1; 0 0 0 1 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0; 0 0 0 0 -3 3 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0 0 0 0 -3 3 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0 0 0 0 -4 4 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0 0 0 0 2 -2 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; 0 0 0 0 3 -3 0  0 0 0 0 0 0 0 0 0 0 0 0 1 0]

#Part B: Determine if the Urea Cycle is elementary balanced

A = [4 7 1 4 0 0; 6 13 3 3 0 0; 10 18 4 6 0 0; 4 4 0 4 0 0; 6 14 4 2 0 0; 1 4 2 1 0 0; 5 12 2 2 0 0; 1 4 1 5 1 0; 10 16 5 13 3 0; 10 14 5 7 1 0; 0 4 0 7 2 0; 0 2 0 1 0 0; 0 3 0 4 1 0; 21 30 7 17 3 0; 0 1 0 0 0 0; 0 0 0 2 0 0; 0 0 1 1 0 0; 21 29 7 17 3 0]
#A is the atom matrix
Check_Balance1 = transpose(A)*S #If first six columns are zero then it is elementally balanced

#Part C: Calculate Maximum rate of Urea Production
include(“Flux.jl”)
Kcat1 = 203 #sec^-1
Kcat2 = 34.5 
Kcat3 = 13.7
Kcat4 = 88.1
Kcat5 = 249

E = 0.01e-3 #mmol/gDW
E_new = E * ((2.3e-9)/((3.7e-12)*0.7)) #mM

#Km values found using BRENDA
Km_v1_Asp = 0.017 #mM
Km_v1_ATP = 0.041
Km_v3_Arg = 2.5
Km_v5_Arg = 2.5
Km_v5_NADPH = 0.00525

#Metabolic Concentrations from Park et al. Paper
ATP_v1 = 4.67 #mM
Asp_v1 = 14.9
Arg_v3 = 2.55e-1
NADPH_v5 = 6.54e-2
Arg_v5 = 2.55e-1

#Saturation calculated with the general form of x/(x + Km)
Saturation_1 = (Asp_v1 / (Asp_v1 + Km_v1_Asp))*(ATP_v1 / (ATP_v1 + Km_v1_ATP))
Saturation_3 = (Arg_v3 / (Arg_v3 + Km_v3_Arg))
Saturation_5 = (Arg_v5 / (Arg_v5 + Km_v5_Arg))*(NADPH_v5 / (NADPH_v5 + Km_v5_NADPH))

#Calculate Upper Bounds (UB)
UB1 = Kcat1*Saturation_1*E_new #mmol/gDW-s
UB2 = Kcat2*E_new
UB3 = Kcat3*Saturation_3*E_new
UB4 = Kcat4*E_new
UB5 = Kcat5*Saturation_5*E_new
UBB = (10*((2.3e-9)/((3.7e-12)*0.7)))/3600 #mM/s

Flux_Array = [0 UB1; 0 UB2; 0 UB3; 0 UB4; 0 UB5; 0 UB5; 0 UBB; 0 UBB; 0 UBB; 0 UBB; 0 UBB; 0 UBB; 0 UBB; 0 UBB; 0 UBB; 0 UBB; 0 UBB; 0 UBB; 0 UBB; 0 UBB; 0 UBB] 

Objective_Array = [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

Species_Array = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]

(Objective_value, flux_array, uptake_array, exit_flag) = calculate_optimal_flux_distribution(S, Flux_Array, Species_Bound, Objective_Array)

Maxrate = -1*objective_value
