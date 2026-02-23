#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import bokeh.io
import bokeh.plotting
import scipy 
import colorcet
import panel as pn
import math
from bokeh.io import export_png
import holoviews as hv
from bokeh.io import export_svg
from sympy.solvers import solve

bokeh.io.output_notebook()
# Set up color palette for this notebook
colors = colorcet.b_glasbey_category10

pn.extension()
Rim4_tau_slider = pn.widgets.FloatSlider(
    name="Rim4 τ", start=0, end=10, step=0.01, value=0.2
)
Rim4_T_slider = pn.widgets.FloatSlider(
    name="Rim4 T", start=20, end=100, step=1, value=20
)

@pn.depends(
    Rim4_tau_slider.param.value,
    Rim4_T_slider.param.value,
   
)

def Rim4_plot(
    Rim4_tau=0.2,
    Rim4_T=20,
    t_max=140,
):
    args = (
        Rim4_tau,
        Rim4_T,
    )
  
    t = np.linspace(0, t_max, 400)
    Rim4 = (1-np.tanh(Rim4_tau*(t-Rim4_T)))/2
    Rim4Slow = (1-np.tanh(Rim4_tau*(t-28)))/2


    # Set up plot
    p = bokeh.plotting.figure(
        x_axis_label="Time (minutes)",
        y_axis_label="mRNA-Rim4 Complex Concentration (a.u.)",
        x_range=[0, t_max],
    )

      # Populate glyphs
    p.line(t, Rim4, line_width=2, color=colors[1], legend_label="Wildtype")
    p.line(t, Rim4Slow, line_width=2, color=colors[69], legend_label="Slow Degradation")

    

    # Place the legend
    p.legend.location = "top_right"
    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"
    p.yaxis.major_label_text_font_size = '12pt'
    p.xaxis.major_label_text_font_size = '11pt'

    return p

widgets = pn.Column(
    pn.Spacer(height=10),
    Rim4_tau_slider,
    Rim4_T_slider,
    width=150,
)

left_column = pn.Column(
    pn.Row(pn.Spacer(width=30)), Rim4_plot,
)

# Final layout
pn.Row(left_column, pn.Spacer(width=20), widgets)



# In[5]:


#--------Rate constants---------
#(In terms of protein numbers: )

   

V = 10**(-11); 
NA = 6.023*10**(23);
vna = V*NA*(100*10**(-9))*10**(-3); 
ivna = 1/vna;

#(Clb4)

kClb4s = 0.0008;
kClb4sp = 0.009;
kClb4d = 0.02;
kClb4dp = 0.2;
kClb4dpp = 0.02;

#(Clb1)

kClb1s = 0.04;
kClb1sp = 0.024;
kClb1d = 0.1;
kClb1dp = 0.5;
kClb1dpp = 0.02;

#(Clb3)

kClb3s = 0.002;
kClb3sp = 0.000114;
kClb3d = 0.34;
kClb3dp = 0.2;
kClb3dpp = 0.02;
AlphaClb3 = 0.4;


#(Cdc20)

kCdc20s = 0.1;
kCdc20sp = 0.001;
kCdc20d = 0.05;
kCdc20dp = 0.02; 
kClb1Cdc20d = 0.8;
kClb4Cdc20d = 0.2;
kClb3Cdc20d = 0.4;


#(New model in progress: Includes APC, Cdc20 phosophorylation and binding )

kApcCdc20a = 0.25;
kApcCdc20d = 0.32;
JApcCdc20a = 0.1;
JApcCdc20d = 0.1;
kApcClb = 0.09;
JApcClb = 0.001;
kApcp = 0.072; 
JApcp = 0.01;
APCtot = 5*vna;

#(Cdc5)

kCdc5a = 0.1;
kCdc5i = 2;
kCdc5d = 0.04;
kCdc5dp = 0.003;
kCdc5dpp = 0.002;
kCdc5ap = 1.2;
kCdc5app = 2;
kCdc5s = 0.004; kCdc5sp = 0.027;

#(Ndt80)

kNdt80s = 0.01;
kNdt80sp = 2;
kNdt80d = 0.093;
kNdt80dp = 0.013; 
n_Ndt80 = 1.62;
K_Ndt80  = 6338;

#(Ama1)

kAma1i = 0.005;
kAma1ip = 0.5;
JAma1 = 0.1;
kAma1a = 0.1;
kAma1s = 0.01;
kAma1dp = 0.02;
AlphakAma1s = 0.4;


def Rim4(t,Rim4_tau, Rim4_T):
    R4 = (1-np.tanh(Rim4_tau*(t-Rim4_T)))/2
    return R4
  
def Meiosis_ODE(Odes, t, Rim4_tau, Rim4_T):
    
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N = Odes


    # Compute dOdes/dt
    dClb1_dt = 1.0*(kClb1s*vna + kClb1sp * Ndt80N - (kClb1d +  kClb1dp*ivna * Ama1N + kClb1dpp*ivna * (Ama1N + Ama1PN) + kClb1Cdc20d*ivna*APCpCdc20N) * Clb1N)
    dClb4_dt = 1.0*(kClb4s*vna + kClb4sp * Ndt80N - (kClb4d +  kClb4dp*ivna * Ama1N + kClb4dpp*ivna * (Ama1N + Ama1PN) + kClb4Cdc20d*ivna*APCpCdc20N) * Clb4N)
    dClb3_dt = 1.0*(kClb3s*vna - (kClb3d  + (kClb3dp*ivna * Ama1N + kClb3dpp*ivna * (Ama1N + Ama1PN)) + kClb3Cdc20d*ivna*APCpCdc20N) * Clb3N + AlphaClb3*(vna - Rim4(t,Rim4_tau, Rim4_T)*vna)*math.exp(-(t - Rim4_T)/25) + kClb3sp*(1 - Rim4(t,Rim4_tau, Rim4_T))* Ndt80N) 
    dNdt80_dt = 1.0*(kNdt80s*vna + (kNdt80sp*vna*(Ndt80N)**(n_Ndt80))/(K_Ndt80**(n_Ndt80) + Ndt80N**(n_Ndt80)) - kNdt80d * Ndt80N - kNdt80dp*ivna*Ama1N* Ndt80N) 
    dCdc20T_dt = 1.0*(kCdc20s*vna + kCdc20sp * Ndt80N - kCdc20d*Cdc20TN - kCdc20dp* APCpCdc20N)
    dAPCp_dt = kApcClb*Clb4N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) + kApcClb*Clb1N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) - kApcp* vna*APCpN/(JApcp*vna + APCpN) + kApcClb*Clb3N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N))
    dAPCpCdc20_dt = kApcCdc20a*APCpN*(Cdc20TN - APCpCdc20N)/(JApcCdc20a*vna + Cdc20TN - APCpCdc20N) - kApcCdc20d*vna*APCpCdc20N/(JApcCdc20d*vna + APCpCdc20N) 
    dCdc5A_dt = 1.0*((kCdc5a + kCdc5ap *ivna* Clb1N + kCdc5ap *ivna* Clb4N)*(Cdc5TN - Cdc5AN) - kCdc5i*Cdc5AN - (kCdc5d + kCdc5dp* ivna*Ama1N + kCdc5dpp* ivna*(Ama1N + Ama1PN) ) * Cdc5AN) + kCdc5app *ivna* Clb3N*(Cdc5TN - Cdc5AN)
    dCdc5T_dt = 1.0*(kCdc5s*vna + kCdc5sp * Ndt80N - (kCdc5d + kCdc5dp * ivna*Ama1N + kCdc5dpp * ivna*(Ama1N + Ama1PN)) * Cdc5TN)
    dAma1P_dt = 1.0*((kAma1i + kAma1ip*ivna*Clb1N + kAma1ip*ivna*Clb4N + kAma1ip*ivna*Clb3N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - kAma1dp*Ama1PN)
    dAma1_dt = 1.0*(kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - (kAma1i + kAma1ip*ivna*Clb1N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb4N * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb3N * Ama1N*vna/(JAma1*vna + Ama1N) + kAma1s*vna - kAma1dp*Ama1N + AlphakAma1s*(1*vna - Rim4(t,Rim4_tau, Rim4_T)*vna)*math.exp(-(t - Rim4_T)/25))

    # Return the result as a NumPy array
    return np.array([dClb1_dt,  dClb4_dt, dClb3_dt, dNdt80_dt, dCdc20T_dt, dAPCp_dt, dAPCpCdc20_dt, dCdc5A_dt,  dCdc5T_dt, dAma1P_dt, dAma1_dt])

# Initial conditions
Odes_0 = np.array([1.125* vna,0.12* vna,0.0,5 *vna,2*vna,0.0,0.1*vna,0.25 *vna,0.75*vna,0.0,0.0])

#Iteractive parameters definition 
pn.extension()
Rim4_tau_slider = pn.widgets.FloatSlider(
    name="Rim4 τ", start=0, end=10, step=0.01, value=0.2
)
Rim4_T_slider = pn.widgets.FloatSlider(
    name="Rim4 T", start=1, end=200, step=1, value=20
)
t_max_slider = pn.widgets.FloatSlider(
    name="t max",start=1, end=300, step=1, value=140
)
@pn.depends(
    Rim4_tau_slider.param.value,
    Rim4_T_slider.param.value,
    t_max_slider.param.value,
   
)

#Plotting function definition 
def WildType_plot(
    Rim4_tau=0.2,
    Rim4_T=20,
    t_max=300,
):
    args = (
        Rim4_tau,
        Rim4_T
    )

      # Integrate ODES
    t = np.linspace(0, t_max, 400)
    Odes = scipy.integrate.odeint(Meiosis_ODE, Odes_0, t, args=args)
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N  = Odes.transpose()
  
    # Set up plot
    p = bokeh.plotting.figure(
        frame_width=425,
        frame_height=350,
        x_axis_label="Time (minutes)",
        y_axis_label="Protein Concentration (a.u.)",
        x_range=[0, t_max],
        
    )

      # Plot of protein concentrations in arbitrary units 
    p.line(t, Clb1N/vna, line_width=2, color=colors[37], legend_label="Clb1")
    p.line(t, Clb4N/vna, line_width=2, color=colors[20], legend_label="Clb4")
    p.line(t, Clb3N/vna, line_width=2, color=colors[72], legend_label="Clb3")
    p.line(t, Ndt80N/vna, line_width=2, color=colors[3], legend_label="Ndt80N")
    p.line(t, Cdc20TN/vna, line_width=2, color=colors[4], legend_label="Cdc20TN")
    p.line(t, APCpN/vna, line_width=2, color=colors[87], legend_label="APCpN")
    p.line(t, APCpCdc20N/vna, line_width=2, color=colors[9], legend_label="APCpCdc20")
    p.line(t, Cdc5AN/vna, line_width=2, color=colors[70], legend_label="Cdc5A")
    p.line(t, Cdc5TN/vna, line_width=2, color=colors[8], legend_label="Cdc5TN")
    p.line(t, Ama1PN/vna, line_width=2, color=colors[14], legend_label="Ama1PN")
    p.line(t, Ama1N/vna , line_width=2, color=colors[10], legend_label="Ama1N")
    p.line(t, (Clb1N+Clb3N+Clb4N)/vna, line_width=3, color=colors[11], legend_label="Clb Total")
    p.line(t,0.47,line_width=3,color=colors[11],line_dash='4 4')

    

    # Place the legend
    p.legend.location = "top_right"
    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"
    p.yaxis.major_label_text_font_size = '12pt'
    p.xaxis.major_label_text_font_size = '11pt'
    
    return p
    
widgets = pn.Column(
    pn.Spacer(height=10),
    Rim4_tau_slider,
    Rim4_T_slider,
    t_max_slider,
    width=150,
)

left_column = pn.Column(
    pn.Row(pn.Spacer(width=20), WildType_plot,
))

# Final layout
pn.Row(left_column, pn.Spacer(width=20),widgets)





# In[9]:


#----Ama1 mutant----

def Ama1Mut_ODE(Odes, t, Rim4_tau, Rim4_T):
    
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N = Odes

    # Compute dOdes/dt
    dClb1_dt = 1.0*(kClb1s*vna + kClb1sp * Ndt80N - (kClb1d +  kClb1dp*ivna * Ama1N + kClb1dpp*ivna * (Ama1N + Ama1PN) + kClb1Cdc20d*ivna*APCpCdc20N) * Clb1N)
    dClb4_dt = 1.0*(kClb4s*vna + kClb4sp * Ndt80N - (kClb4d +  kClb4dp*ivna * Ama1N + kClb4dpp*ivna * (Ama1N + Ama1PN) + kClb4Cdc20d*ivna*APCpCdc20N) * Clb4N)
    dClb3_dt = 1.0*(kClb3s*vna - (kClb3d  + (kClb3dp*ivna * Ama1N + kClb3dpp*ivna * (Ama1N + Ama1PN)) + kClb3Cdc20d*ivna*APCpCdc20N) * Clb3N + AlphaClb3*(vna - Rim4(t,Rim4_tau, Rim4_T)*vna)*math.exp(-(t - Rim4_T)/25) + kClb3sp*(1 - Rim4(t,Rim4_tau, Rim4_T))* Ndt80N) 
    dNdt80_dt = 1.0*(kNdt80s*vna + (kNdt80sp*vna*(Ndt80N)**(n_Ndt80))/(K_Ndt80**(n_Ndt80) + Ndt80N**(n_Ndt80)) - kNdt80d * Ndt80N - kNdt80dp*ivna*Ama1N* Ndt80N) 
    dCdc20T_dt = 1.0*(kCdc20s*vna + kCdc20sp * Ndt80N - kCdc20d*Cdc20TN - kCdc20dp* APCpCdc20N)
    dAPCp_dt = kApcClb*Clb4N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) + kApcClb*Clb1N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) - kApcp* vna*APCpN/(JApcp*vna + APCpN) + kApcClb*Clb3N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N))
    dAPCpCdc20_dt = kApcCdc20a*APCpN*(Cdc20TN - APCpCdc20N)/(JApcCdc20a*vna + Cdc20TN - APCpCdc20N) - kApcCdc20d*vna*APCpCdc20N/(JApcCdc20d*vna + APCpCdc20N) 
    dCdc5A_dt = 1.0*((kCdc5a + kCdc5ap *ivna* Clb1N + kCdc5ap *ivna* Clb4N)*(Cdc5TN - Cdc5AN) - kCdc5i*Cdc5AN - (kCdc5d + kCdc5dp* ivna*Ama1N + kCdc5dpp* ivna*(Ama1N + Ama1PN) ) * Cdc5AN) + kCdc5app *ivna* Clb3N*(Cdc5TN - Cdc5AN)
    dCdc5T_dt = 1.0*(kCdc5s*vna + kCdc5sp * Ndt80N - (kCdc5d + kCdc5dp * ivna*Ama1N + kCdc5dpp * ivna*(Ama1N + Ama1PN)) * Cdc5TN)
    dAma1P_dt = 0.0*((kAma1i + kAma1ip*ivna*Clb1N + kAma1ip*ivna*Clb4N + kAma1ip*ivna*Clb3N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - kAma1dp*Ama1PN)
    dAma1_dt = 0.0*(kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - (kAma1i + kAma1ip*ivna*Clb1N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb4N * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb3N * Ama1N*vna/(JAma1*vna + Ama1N) + kAma1s*vna - kAma1dp*Ama1N + AlphakAma1s*(1*vna - Rim4(t,Rim4_tau, Rim4_T)*vna)*math.exp(-(t - Rim4_T)/25))

    # Return the result as a NumPy array
    return np.array([dClb1_dt,  dClb4_dt, dClb3_dt, dNdt80_dt, dCdc20T_dt, dAPCp_dt, dAPCpCdc20_dt, dCdc5A_dt,  dCdc5T_dt, dAma1P_dt, dAma1_dt])

# Initial condition
Odes_0 = np.array([1.125* vna,0.12* vna,0.0,5 *vna,2*vna,0.0,0.1*vna,0.25 *vna,0.75*vna,0.0,0.0])

#Sliders

pn.extension()
Rim4_tau_slider = pn.widgets.FloatSlider(
    name="Rim4 τ", start=0, end=10, step=0.01, value=0.2
)
Rim4_T_slider = pn.widgets.FloatSlider(
    name="Rim4 T", start=1, end=100, step=1, value=20
)
t_max_slider = pn.widgets.FloatSlider(
    name="t max",start=1, end=300, step=1, value=140
)
@pn.depends(
    Rim4_tau_slider.param.value,
    Rim4_T_slider.param.value,
    t_max_slider.param.value,
   
)

def Ama1Mut_plot(
    Rim4_tau=0.2,
    Rim4_T=20,
    t_max=200,
):
    args = (
        Rim4_tau,
        Rim4_T
    )

      # Integrate ODES
    t = np.linspace(0, t_max, 400)
    Odes = scipy.integrate.odeint(Ama1Mut_ODE, Odes_0, t, args=args)
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N  = Odes.transpose()
  
    # Set up plot
    p = bokeh.plotting.figure(
        frame_width=425,
        frame_height=350,
        x_axis_label="Time (minutes)",
        y_axis_label="Protein Concentration (a.u.)",
        x_range=[0, t_max],
    )

      # Populate glyphs
    p.line(t, Clb1N/vna, line_width=2, color=colors[37], legend_label="Clb1")
    p.line(t, Clb4N/vna, line_width=2, color=colors[20], legend_label="Clb4")
    p.line(t, Clb3N/vna, line_width=2, color=colors[72], legend_label="Clb3")
    p.line(t, Ndt80N/vna, line_width=2, color=colors[3], legend_label="Ndt80N")
    p.line(t, Cdc20TN/vna, line_width=2, color=colors[4], legend_label="Cdc20TN")
    p.line(t, APCpN/vna, line_width=2, color=colors[5], legend_label="APCpN")
    p.line(t, APCpCdc20N/vna, line_width=2, color=colors[9], legend_label="APCpCdc20")
    p.line(t, Cdc5AN/vna, line_width=2, color=colors[70], legend_label="Cdc5A")
    p.line(t, Cdc5TN/vna, line_width=2, color=colors[8], legend_label="Cdc5TN")
    p.line(t, Ama1PN/vna, line_width=2, color=colors[9], legend_label="Ama1PN")
    p.line(t, Ama1N/vna, line_width=2, color=colors[10], legend_label="Ama1N")
    p.line(t, (Clb1N+Clb3N+Clb4N)/vna, line_width=2, color=colors[11], legend_label="Clb Total")
    p.line(t,0.47,line_width=3,color=colors[11],line_dash='4 4')

    

     # Place the legend
    p.legend.location = "top_right"
    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"
    p.yaxis.major_label_text_font_size = '12pt'
    p.xaxis.major_label_text_font_size = '11pt'
    


    return p
    
widgets = pn.Column(
    pn.Spacer(height=10),
    Rim4_tau_slider,
    Rim4_T_slider,
    t_max_slider,
    width=150,
)
left_column = pn.Column(
    pn.Row(pn.Spacer(width=10)), Ama1Mut_plot,
)

# Final layout
pn.Row(left_column, pn.Spacer(width=20),widgets)


# In[11]:


#----Clb1 mutant----

def Clb1Mut_ODE(Odes, t, Rim4_tau, Rim4_T):
    
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N = Odes

    # Compute dOdes/dt
    dClb1_dt = 0.0*(kClb1s*vna + kClb1sp * Ndt80N - (kClb1d +  kClb1dp*ivna * Ama1N + kClb1dpp*ivna * (Ama1N + Ama1PN) + kClb1Cdc20d*ivna*APCpCdc20N) * Clb1N)
    dClb4_dt = 1.0*(kClb4s*vna + kClb4sp * Ndt80N - (kClb4d +  kClb4dp*ivna * Ama1N + kClb4dpp*ivna * (Ama1N + Ama1PN) + kClb4Cdc20d*ivna*APCpCdc20N) * Clb4N)
    dClb3_dt = 1.0*(kClb3s*vna - (kClb3d  + (kClb3dp*ivna * Ama1N + kClb3dpp*ivna * (Ama1N + Ama1PN)) + kClb3Cdc20d*ivna*APCpCdc20N) * Clb3N + AlphaClb3*(vna - Rim4(t,Rim4_tau, Rim4_T)*vna)*math.exp(-(t - Rim4_T)/25) + kClb3sp*(1 - Rim4(t,Rim4_tau, Rim4_T))* Ndt80N) 
    dNdt80_dt = 1.0*(kNdt80s*vna + (kNdt80sp*vna*(Ndt80N)**(n_Ndt80))/(K_Ndt80**(n_Ndt80) + Ndt80N**(n_Ndt80)) - kNdt80d * Ndt80N - kNdt80dp*ivna*Ama1N* Ndt80N) 
    dCdc20T_dt = 1.0*(kCdc20s*vna + kCdc20sp * Ndt80N - kCdc20d*Cdc20TN - kCdc20dp* APCpCdc20N)
    dAPCp_dt = kApcClb*Clb4N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) + kApcClb*Clb1N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) - kApcp* vna*APCpN/(JApcp*vna + APCpN) + kApcClb*Clb3N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N))
    dAPCpCdc20_dt = kApcCdc20a*APCpN*(Cdc20TN - APCpCdc20N)/(JApcCdc20a*vna + Cdc20TN - APCpCdc20N) - kApcCdc20d*vna*APCpCdc20N/(JApcCdc20d*vna + APCpCdc20N) 
    dCdc5A_dt = 1.0*((kCdc5a + kCdc5ap *ivna* Clb1N + kCdc5ap *ivna* Clb4N)*(Cdc5TN - Cdc5AN) - kCdc5i*Cdc5AN - (kCdc5d + kCdc5dp* ivna*Ama1N + kCdc5dpp* ivna*(Ama1N + Ama1PN) ) * Cdc5AN) + kCdc5app *ivna* Clb3N*(Cdc5TN - Cdc5AN)
    dCdc5T_dt = 1.0*(kCdc5s*vna + kCdc5sp * Ndt80N - (kCdc5d + kCdc5dp * ivna*Ama1N + kCdc5dpp * ivna*(Ama1N + Ama1PN)) * Cdc5TN)
    dAma1P_dt = 1.0*((kAma1i + kAma1ip*ivna*Clb1N + kAma1ip*ivna*Clb4N + kAma1ip*ivna*Clb3N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - kAma1dp*Ama1PN)
    dAma1_dt = 1.0*(kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - (kAma1i + kAma1ip*ivna*Clb1N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb4N * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb3N * Ama1N*vna/(JAma1*vna + Ama1N) + kAma1s*vna - kAma1dp*Ama1N + AlphakAma1s*(1*vna - Rim4(t,Rim4_tau, Rim4_T)*vna)*math.exp(-(t - Rim4_T)/25))

    # Return the result as a NumPy array
    return np.array([dClb1_dt,  dClb4_dt, dClb3_dt, dNdt80_dt, dCdc20T_dt, dAPCp_dt, dAPCpCdc20_dt, dCdc5A_dt,  dCdc5T_dt, dAma1P_dt, dAma1_dt])

# Initial condition
Odes_0 = np.array([0,0.12* vna,0.0,5 *vna,2*vna,0.0,0.1*vna,0.25 *vna,0.75*vna,0.0,0.0])

#Sliders

pn.extension()
Rim4_tau_slider = pn.widgets.FloatSlider(
    name="Rim4 τ", start=0, end=10, step=0.01, value=0.2
)
Rim4_T_slider = pn.widgets.FloatSlider(
    name="Rim4 T", start=1, end=100, step=1, value=20
)
t_max_slider = pn.widgets.FloatSlider(
    name="t max",start=1, end=300, step=1, value=140
)
@pn.depends(
    Rim4_tau_slider.param.value,
    Rim4_T_slider.param.value,
    t_max_slider.param.value,
   
)

def Clb1Mut_plot(
    Rim4_tau=0.2,
    Rim4_T=20,
    t_max=140,
):
    args = (
        Rim4_tau,
        Rim4_T
    )

      # Integrate ODES
    t = np.linspace(0, t_max, 400)
    Odes = scipy.integrate.odeint(Clb1Mut_ODE, Odes_0, t, args=args)
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N  = Odes.transpose()
  
    # Set up plot
    p = bokeh.plotting.figure(
        frame_width=425,
        frame_height=350,
        x_axis_label="Time (minutes)",
        y_axis_label="Protein Concentration (a.u.)",
        x_range=[0, t_max],
    )

      # Populate glyphs
    p.line(t, Clb1N/vna, line_width=2, color=colors[37], legend_label="Clb1")
    p.line(t, Clb4N/vna, line_width=2, color=colors[20], legend_label="Clb4")
    p.line(t, Clb3N/vna, line_width=2, color=colors[72], legend_label="Clb3")
    p.line(t, Ndt80N/vna, line_width=2, color=colors[3], legend_label="Ndt80N")
    p.line(t, Cdc20TN/vna, line_width=2, color=colors[4], legend_label="Cdc20TN")
    p.line(t, APCpN/vna, line_width=2, color=colors[5], legend_label="APCpN")
    p.line(t, APCpCdc20N/vna, line_width=2, color=colors[9], legend_label="APCpCdc20")
    p.line(t, Cdc5AN/vna, line_width=2, color=colors[70], legend_label="Cdc5A")
    p.line(t, Cdc5TN/vna, line_width=2, color=colors[8], legend_label="Cdc5TN")
    p.line(t, Ama1PN/vna, line_width=2, color=colors[9], legend_label="Ama1PN")
    p.line(t, Ama1N/vna , line_width=2, color=colors[10], legend_label="Ama1N")
    p.line(t, (Clb1N+Clb3N+Clb4N)/vna, line_width=2, color=colors[11], legend_label="Clb Total")
    p.line(t,0.47,line_width=3,color=colors[11],line_dash='4 4')
    

    

     # Place the legend
    p.legend.location = "top_right"
    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"
    p.yaxis.major_label_text_font_size = '12pt'
    p.xaxis.major_label_text_font_size = '11pt'


    return p
    
widgets = pn.Column(
    pn.Spacer(height=10),
    Rim4_tau_slider,
    Rim4_T_slider,
    t_max_slider,
    width=150,
)

left_column = pn.Column(
    pn.Row(pn.Spacer(width=10)), Clb1Mut_plot,
)

# Final layout
pn.Row(left_column, pn.Spacer(width=20),widgets)


# In[13]:


#----Clb3 mutant----

def Clb3Mut_ODE(Odes, t, Rim4_tau, Rim4_T):
    
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N = Odes

    # Compute dOdes/dt
    dClb1_dt = 1.0*(kClb1s*vna + kClb1sp * Ndt80N - (kClb1d +  kClb1dp*ivna * Ama1N + kClb1dpp*ivna * (Ama1N + Ama1PN) + kClb1Cdc20d*ivna*APCpCdc20N) * Clb1N)
    dClb4_dt = 1.0*(kClb4s*vna + kClb4sp * Ndt80N - (kClb4d +  kClb4dp*ivna * Ama1N + kClb4dpp*ivna * (Ama1N + Ama1PN) + kClb4Cdc20d*ivna*APCpCdc20N) * Clb4N)
    dClb3_dt = 0.0*(kClb3s*vna - (kClb3d  + (kClb3dp*ivna * Ama1N + kClb3dpp*ivna * (Ama1N + Ama1PN)) + kClb3Cdc20d*ivna*APCpCdc20N) * Clb3N + AlphaClb3*(vna - Rim4(t,Rim4_tau, Rim4_T)*vna)*math.exp(-(t - Rim4_T)/25) + kClb3sp*(1 - Rim4(t,Rim4_tau, Rim4_T))* Ndt80N) 
    dNdt80_dt = 1.0*(kNdt80s*vna + (kNdt80sp*vna*(Ndt80N)**(n_Ndt80))/(K_Ndt80**(n_Ndt80) + Ndt80N**(n_Ndt80)) - kNdt80d * Ndt80N - kNdt80dp*ivna*Ama1N* Ndt80N)  
    dCdc20T_dt = 1.0*(kCdc20s*vna + kCdc20sp * Ndt80N - kCdc20d*Cdc20TN - kCdc20dp* APCpCdc20N)
    dAPCp_dt = kApcClb*Clb4N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) + kApcClb*Clb1N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) - kApcp* vna*APCpN/(JApcp*vna + APCpN) + kApcClb*Clb3N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N))
    dAPCpCdc20_dt = kApcCdc20a*APCpN*(Cdc20TN - APCpCdc20N)/(JApcCdc20a*vna + Cdc20TN - APCpCdc20N) - kApcCdc20d*vna*APCpCdc20N/(JApcCdc20d*vna + APCpCdc20N) 
    dCdc5A_dt = 1.0*((kCdc5a + kCdc5ap *ivna* Clb1N + kCdc5ap *ivna* Clb4N)*(Cdc5TN - Cdc5AN) - kCdc5i*Cdc5AN - (kCdc5d + kCdc5dp* ivna*Ama1N + kCdc5dpp* ivna*(Ama1N + Ama1PN) ) * Cdc5AN) + kCdc5app *ivna* Clb3N*(Cdc5TN - Cdc5AN)
    dCdc5T_dt = 1.0*(kCdc5s*vna + kCdc5sp * Ndt80N - (kCdc5d + kCdc5dp * ivna*Ama1N + kCdc5dpp * ivna*(Ama1N + Ama1PN)) * Cdc5TN)
    dAma1P_dt = 1.0*((kAma1i + kAma1ip*ivna*Clb1N + kAma1ip*ivna*Clb4N + kAma1ip*ivna*Clb3N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - kAma1dp*Ama1PN)
    dAma1_dt = 1.0*(kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - (kAma1i + kAma1ip*ivna*Clb1N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb4N * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb3N * Ama1N*vna/(JAma1*vna + Ama1N) + kAma1s*vna - kAma1dp*Ama1N + AlphakAma1s*(1*vna - Rim4(t,Rim4_tau, Rim4_T)*vna)*math.exp(-(t - Rim4_T)/25))

    # Return the result as a NumPy array
    return np.array([dClb1_dt,  dClb4_dt, dClb3_dt, dNdt80_dt, dCdc20T_dt, dAPCp_dt, dAPCpCdc20_dt, dCdc5A_dt,  dCdc5T_dt, dAma1P_dt, dAma1_dt])

# Initial condition
Odes_0 = np.array([1.125* vna,0.12* vna,0.0,5 *vna,2*vna,0.0,0.1*vna,0.25 *vna,0.75*vna,0.0,0.0])

#Sliders

pn.extension()
Rim4_tau_slider = pn.widgets.FloatSlider(
    name="Rim4 τ", start=0, end=10, step=0.01, value=0.2
)
Rim4_T_slider = pn.widgets.FloatSlider(
    name="Rim4 T", start=1, end=100, step=1, value=20
)
t_max_slider = pn.widgets.FloatSlider(
    name="t max",start=1, end=300, step=1, value=140
)
@pn.depends(
    Rim4_tau_slider.param.value,
    Rim4_T_slider.param.value,
    t_max_slider.param.value,
   
)
def Clb3Mut_plot(
    Rim4_tau=0.2,
    Rim4_T=20,
    t_max=200,
):
    args = (
        Rim4_tau,
        Rim4_T
    )

      # Integrate ODES
    t = np.linspace(0, t_max, 400)
    Odes = scipy.integrate.odeint(Clb3Mut_ODE, Odes_0, t, args=args)
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N  = Odes.transpose()
  
    # Set up plot
    p = bokeh.plotting.figure(
        frame_width=425,
        frame_height=350,
        x_axis_label="Time (minutes)",
        y_axis_label="Protein Concentration (a.u.)",
        x_range=[0, t_max],
    )

      # Populate glyphs
    p.line(t, Clb1N/vna, line_width=2, color=colors[37], legend_label="Clb1")
    p.line(t, Clb4N/vna, line_width=2, color=colors[20], legend_label="Clb4")
    p.line(t, Clb3N/vna, line_width=2, color=colors[72], legend_label="Clb3")
    p.line(t, Ndt80N/vna, line_width=2, color=colors[3], legend_label="Ndt80N")
    p.line(t, Cdc20TN/vna, line_width=2, color=colors[4], legend_label="Cdc20TN")
    p.line(t, APCpN/vna, line_width=2, color=colors[5], legend_label="APCpN")
    p.line(t, APCpCdc20N/vna, line_width=2, color=colors[9], legend_label="APCpCdc20")
    p.line(t, Cdc5AN/vna, line_width=2, color=colors[70], legend_label="Cdc5A")
    p.line(t, Cdc5TN/vna, line_width=2, color=colors[8], legend_label="Cdc5TN")
    p.line(t, Ama1PN/vna, line_width=2, color=colors[9], legend_label="Ama1PN")
    p.line(t, Ama1N/vna, line_width=2, color=colors[10], legend_label="Ama1N")
    p.line(t, (Clb1N+Clb3N+Clb4N)/vna, line_width=2, color=colors[11], legend_label="Clb Total")
    p.line(t,0.47,line_width=3,color=colors[11],line_dash='4 4')

    

     # Place the legend
    p.legend.location = "top_right"
    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"
    p.yaxis.major_label_text_font_size = '12pt'
    p.xaxis.major_label_text_font_size = '11pt'


    return p
    
widgets = pn.Column(
    pn.Spacer(height=10),
    Rim4_tau_slider,
    Rim4_T_slider,
    t_max_slider,
    width=150,
)

left_column = pn.Column(
    pn.Row(pn.Spacer(width=10)), Clb3Mut_plot,
)

# Final layout
pn.Row(left_column, pn.Spacer(width=20),widgets)


# In[15]:


#----Clb4 mutant----

def Clb4Mut_ODE(Odes, t, Rim4_tau, Rim4_T):
    
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N = Odes

   # Compute dOdes/dt
    dClb1_dt = 1.0*(kClb1s*vna + kClb1sp * Ndt80N - (kClb1d +  kClb1dp*ivna * Ama1N + kClb1dpp*ivna * (Ama1N + Ama1PN) + kClb1Cdc20d*ivna*APCpCdc20N) * Clb1N)
    dClb4_dt = 0.0*(kClb4s*vna + kClb4sp * Ndt80N - (kClb4d +  kClb4dp*ivna * Ama1N + kClb4dpp*ivna * (Ama1N + Ama1PN) + kClb4Cdc20d*ivna*APCpCdc20N) * Clb4N)
    dClb3_dt = 1.0*(kClb3s*vna - (kClb3d  + (kClb3dp*ivna * Ama1N + kClb3dpp*ivna * (Ama1N + Ama1PN)) + kClb3Cdc20d*ivna*APCpCdc20N) * Clb3N + AlphaClb3*(vna - Rim4(t,Rim4_tau, Rim4_T)*vna)*math.exp(-(t - Rim4_T)/25) + kClb3sp*(1 - Rim4(t,Rim4_tau, Rim4_T))* Ndt80N) 
    dNdt80_dt = 1.0*(kNdt80s*vna + (kNdt80sp*vna*(Ndt80N)**(n_Ndt80))/(K_Ndt80**(n_Ndt80) + Ndt80N**(n_Ndt80)) - kNdt80d * Ndt80N - kNdt80dp*ivna*Ama1N* Ndt80N) 
    dCdc20T_dt = 1.0*(kCdc20s*vna + kCdc20sp * Ndt80N - kCdc20d*Cdc20TN - kCdc20dp* APCpCdc20N)
    dAPCp_dt = kApcClb*Clb4N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) + kApcClb*Clb1N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) - kApcp* vna*APCpN/(JApcp*vna + APCpN) + kApcClb*Clb3N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N))
    dAPCpCdc20_dt = kApcCdc20a*APCpN*(Cdc20TN - APCpCdc20N)/(JApcCdc20a*vna + Cdc20TN - APCpCdc20N) - kApcCdc20d*vna*APCpCdc20N/(JApcCdc20d*vna + APCpCdc20N) 
    dCdc5A_dt = 1.0*((kCdc5a + kCdc5ap *ivna* Clb1N + kCdc5ap *ivna* Clb4N)*(Cdc5TN - Cdc5AN) - kCdc5i*Cdc5AN - (kCdc5d + kCdc5dp* ivna*Ama1N + kCdc5dpp* ivna*(Ama1N + Ama1PN) ) * Cdc5AN) + kCdc5app *ivna* Clb3N*(Cdc5TN - Cdc5AN)
    dCdc5T_dt = 1.0*(kCdc5s*vna + kCdc5sp * Ndt80N - (kCdc5d + kCdc5dp * ivna*Ama1N + kCdc5dpp * ivna*(Ama1N + Ama1PN)) * Cdc5TN)
    dAma1P_dt = 1.0*((kAma1i + kAma1ip*ivna*Clb1N + kAma1ip*ivna*Clb4N + kAma1ip*ivna*Clb3N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - kAma1dp*Ama1PN)
    dAma1_dt = 1.0*(kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - (kAma1i + kAma1ip*ivna*Clb1N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb4N * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb3N * Ama1N*vna/(JAma1*vna + Ama1N) + kAma1s*vna - kAma1dp*Ama1N + AlphakAma1s*(1*vna - Rim4(t,Rim4_tau, Rim4_T)*vna)*math.exp(-(t - Rim4_T)/25))
    # Return the result as a NumPy array
    return np.array([dClb1_dt,  dClb4_dt, dClb3_dt, dNdt80_dt, dCdc20T_dt, dAPCp_dt, dAPCpCdc20_dt, dCdc5A_dt,  dCdc5T_dt, dAma1P_dt, dAma1_dt])


# Initial condition
Odes_0 = np.array([1.125* vna,0.0,0.0,5 *vna,2*vna,0.0,0.1*vna,0.25 *vna,0.75*vna,0.0,0.0])
#Sliders

pn.extension()
Rim4_tau_slider = pn.widgets.FloatSlider(
    name="Rim4 τ", start=0, end=10, step=0.01, value=0.2
)
Rim4_T_slider = pn.widgets.FloatSlider(
    name="Rim4 T", start=1, end=100, step=1, value=20
)
t_max_slider = pn.widgets.FloatSlider(
    name="t max",start=1, end=300, step=1, value=140
)
@pn.depends(
    Rim4_tau_slider.param.value,
    Rim4_T_slider.param.value,
    t_max_slider.param.value,
   
)
def Clb4Mut_plot(
    Rim4_tau=0.2,
    Rim4_T=20,
    t_max=200,
):
    args = (
        Rim4_tau,
        Rim4_T
    )

      # Integrate ODES
    t = np.linspace(0, t_max, 400)
    Odes = scipy.integrate.odeint(Clb4Mut_ODE, Odes_0, t, args=args)
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N  = Odes.transpose()
  
    # Set up plot
    p = bokeh.plotting.figure(
        frame_width=425,
        frame_height=350,
        x_axis_label="Time (minutes)",
        y_axis_label="Protein Concentration (a.u.)",
        x_range=[0, t_max],
    )

      # Populate glyphs
    p.line(t, Clb1N/vna, line_width=2, color=colors[37], legend_label="Clb1")
    p.line(t, Clb4N/vna, line_width=2, color=colors[20], legend_label="Clb4")
    p.line(t, Clb3N/vna, line_width=2, color=colors[72], legend_label="Clb3")
    p.line(t, Ndt80N/vna, line_width=2, color=colors[3], legend_label="Ndt80N")
    p.line(t, Cdc20TN/vna, line_width=2, color=colors[4], legend_label="Cdc20TN")
    p.line(t, APCpN/vna, line_width=2, color=colors[5], legend_label="APCpN")
    p.line(t, APCpCdc20N/vna, line_width=2, color=colors[9], legend_label="APCpCdc20")
    p.line(t, Cdc5AN/vna, line_width=2, color=colors[70], legend_label="Cdc5A")
    p.line(t, Cdc5TN/vna, line_width=2, color=colors[8], legend_label="Cdc5TN")
    p.line(t, Ama1PN/vna, line_width=2, color=colors[9], legend_label="Ama1PN")
    p.line(t, Ama1N/vna, line_width=2, color=colors[10], legend_label="Ama1N")
    p.line(t, (Clb1N+Clb3N+Clb4N)/vna, line_width=2, color=colors[11], legend_label="Clb Total")
    p.line(t,0.47,line_width=3,color=colors[11],line_dash='4 4')

    

     # Place the legend
    p.legend.location = "top_right"
    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"
    p.yaxis.major_label_text_font_size = '12pt'
    p.xaxis.major_label_text_font_size = '11pt'


    return p
    
widgets = pn.Column(
    pn.Spacer(height=10),
    Rim4_tau_slider,
    Rim4_T_slider,
    t_max_slider,
    width=150,
)

left_column = pn.Column(
    pn.Row(pn.Spacer(width=10)), Clb4Mut_plot,
)

# Final layout
pn.Row(left_column, pn.Spacer(width=20),widgets)


# In[19]:


#----No Rim4 degradation------
  
def NoRim4Deg_ODE(Odes, t,Rim4_tau,Rim4_T):
    
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N = Odes


    # Compute dOdes/dt
    dClb1_dt = 1.0*(kClb1s*vna + kClb1sp * Ndt80N - (kClb1d +  kClb1dp*ivna * Ama1N + kClb1dpp*ivna * (Ama1N + Ama1PN) + kClb1Cdc20d*ivna*APCpCdc20N) * Clb1N)
    dClb4_dt = 1.0*(kClb4s*vna + kClb4sp * Ndt80N - (kClb4d +  kClb4dp*ivna * Ama1N + kClb4dpp*ivna * (Ama1N + Ama1PN) + kClb4Cdc20d*ivna*APCpCdc20N) * Clb4N)
    dClb3_dt = 1.0*(kClb3s*vna - (kClb3d  + (kClb3dp*ivna * Ama1N + kClb3dpp*ivna * (Ama1N + Ama1PN)) + kClb3Cdc20d*ivna*APCpCdc20N) * Clb3N ) 
    dNdt80_dt = 1.0*(kNdt80s*vna + (kNdt80sp*vna*(Ndt80N)**(n_Ndt80))/(K_Ndt80**(n_Ndt80) + Ndt80N**(n_Ndt80)) - kNdt80d * Ndt80N - kNdt80dp*ivna*Ama1N* Ndt80N) 
    dCdc20T_dt = 1.0*(kCdc20s*vna + kCdc20sp * Ndt80N - kCdc20d*Cdc20TN - kCdc20dp* APCpCdc20N)
    dAPCp_dt = kApcClb*Clb4N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) + kApcClb*Clb1N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N)) - kApcp* vna*APCpN/(JApcp*vna + APCpN) + kApcClb*Clb3N * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N))
    dAPCpCdc20_dt = kApcCdc20a*APCpN*(Cdc20TN - APCpCdc20N)/(JApcCdc20a*vna + Cdc20TN - APCpCdc20N) - kApcCdc20d*vna*APCpCdc20N/(JApcCdc20d*vna + APCpCdc20N) 
    dCdc5A_dt = 1.0*((kCdc5a + kCdc5ap *ivna* Clb1N + kCdc5ap *ivna* Clb4N)*(Cdc5TN - Cdc5AN) - kCdc5i*Cdc5AN - (kCdc5d + kCdc5dp* ivna*Ama1N + kCdc5dpp* ivna*(Ama1N + Ama1PN) ) * Cdc5AN) + kCdc5app *ivna* Clb3N*(Cdc5TN - Cdc5AN)
    dCdc5T_dt = 1.0*(kCdc5s*vna + kCdc5sp * Ndt80N - (kCdc5d + kCdc5dp * ivna*Ama1N + kCdc5dpp * ivna*(Ama1N + Ama1PN)) * Cdc5TN)
    dAma1P_dt = 1.0*((kAma1i + kAma1ip*ivna*Clb1N + kAma1ip*ivna*Clb4N + kAma1ip*ivna*Clb3N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - kAma1dp*Ama1PN)
    dAma1_dt = 1.0*(kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - (kAma1i + kAma1ip*ivna*Clb1N) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb4N * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1ip*ivna*Clb3N * Ama1N*vna/(JAma1*vna + Ama1N) + kAma1s*vna - kAma1dp*Ama1N )

    # Return the result as a NumPy array
    return np.array([dClb1_dt,  dClb4_dt, dClb3_dt, dNdt80_dt, dCdc20T_dt, dAPCp_dt, dAPCpCdc20_dt, dCdc5A_dt,  dCdc5T_dt, dAma1P_dt, dAma1_dt])

# Initial conditions
Odes_0 = np.array([1.125* vna,0.12* vna,0.0,5 *vna,2*vna,0.0,0.1*vna,0.25 *vna,0.75*vna,0.0,0.0])

#Iteractive parameters definition 
pn.extension()
Rim4_tau_slider = pn.widgets.FloatSlider(
    name="Rim4 τ", start=0, end=10, step=0.01, value=0.2
)
Rim4_T_slider = pn.widgets.FloatSlider(
    name="Rim4 T", start=1, end=200, step=1, value=20
)
t_max_slider = pn.widgets.FloatSlider(
    name="t max",start=1, end=300, step=1, value=140
)
@pn.depends(
    Rim4_tau_slider.param.value,
    Rim4_T_slider.param.value,
    t_max_slider.param.value,
   
)

#Plotting function definition 
def NoRim4Deg_plot(
    Rim4_tau=0.2,
    Rim4_T=20,
    t_max=300,
):
    args = (
        Rim4_tau,
        Rim4_T
    )

      # Integrate ODES
    t = np.linspace(0, t_max, 400)
    Odes = scipy.integrate.odeint(NoRim4Deg_ODE, Odes_0, t, args=args)
    Clb1N,  Clb4N, Clb3N, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N  = Odes.transpose()
  
    # Set up plot
    p = bokeh.plotting.figure(
        frame_width=425,
        frame_height=350,
        x_axis_label="Time (minutes)",
        y_axis_label="Protein Concentration (a.u.)",
        x_range=[0, t_max],
        
    )

      # Plot of protein concentrations in arbitrary units 
    p.line(t, Clb1N/vna, line_width=2, color=colors[37], legend_label="Clb1")
    p.line(t, Clb4N/vna, line_width=2, color=colors[20], legend_label="Clb4")
    p.line(t, Clb3N/vna, line_width=2, color=colors[72], legend_label="Clb3")
    p.line(t, Ndt80N/vna, line_width=2, color=colors[3], legend_label="Ndt80N")
    p.line(t, Cdc20TN/vna, line_width=2, color=colors[4], legend_label="Cdc20TN")
    p.line(t, APCpN/vna, line_width=2, color=colors[87], legend_label="APCpN")
    p.line(t, APCpCdc20N/vna, line_width=2, color=colors[9], legend_label="APCpCdc20")
    p.line(t, Cdc5AN/vna, line_width=2, color=colors[70], legend_label="Cdc5A")
    p.line(t, Cdc5TN/vna, line_width=2, color=colors[8], legend_label="Cdc5TN")
    p.line(t, Ama1PN/vna, line_width=2, color=colors[14], legend_label="Ama1PN")
    p.line(t, Ama1N/vna , line_width=2, color=colors[10], legend_label="Ama1N")
    p.line(t, (Clb1N+Clb3N+Clb4N)/vna, line_width=3, color=colors[11], legend_label="Clb Total")
    p.line(t,0.47,line_width=3,color=colors[11],line_dash='4 4')

    

    # Place the legend
    p.legend.location = "top_right"
    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"
    p.yaxis.major_label_text_font_size = '12pt'
    p.xaxis.major_label_text_font_size = '11pt'
    
    return p
    
widgets = pn.Column(
    pn.Spacer(height=10),
    Rim4_tau_slider,
    Rim4_T_slider,
    t_max_slider,
    width=150,
)

left_column = pn.Column(
    pn.Row(pn.Spacer(width=20), NoRim4Deg_plot,
))

# Final layout
pn.Row(left_column, pn.Spacer(width=20),widgets)



# In[27]:


#---- Single ClbT species ------
  
def ClbT_ODE(Odes, t,Rim4_tau,Rim4_T):
    
    ClbTN, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N = Odes


    dClbT_dt = 1.0*(kClb1s*vna + (kClb1sp)*Ndt80N - (kClb1d+(kClb1dp)*ivna*Ama1N + (kClb1dpp)*ivna*(Ama1N + Ama1PN) + (kClb1Cdc20d)*ivna*APCpCdc20N )*ClbTN)
    dNdt80_dt = 1.0*(kNdt80s*vna + (kNdt80sp*vna*(Ndt80N)**(n_Ndt80))/(K_Ndt80**(n_Ndt80) + Ndt80N**(n_Ndt80)) - kNdt80d * Ndt80N - kNdt80dp*ivna*Ama1N* Ndt80N)  
    dCdc20T_dt = 1.0*(kCdc20s*vna + kCdc20sp * Ndt80N - kCdc20d*Cdc20TN - kCdc20dp* APCpCdc20N)
    dAPCp_dt = kApcClb*ClbTN * (APCtot - APCpN - APCpCdc20N)/(JApcClb*vna + (APCtot - APCpN - APCpCdc20N))  - kApcp* vna*APCpN/(JApcp*vna + APCpN) 
    dAPCpCdc20_dt = kApcCdc20a*APCpN*(Cdc20TN - APCpCdc20N)/(JApcCdc20a*vna + Cdc20TN - APCpCdc20N) - kApcCdc20d*vna*APCpCdc20N/(JApcCdc20d*vna + APCpCdc20N) 
    dCdc5A_dt = 1.0*((kCdc5a + kCdc5ap *ivna* ClbTN)*(Cdc5TN - Cdc5AN) - kCdc5i*Cdc5AN - (kCdc5d + kCdc5dp* ivna*Ama1N + kCdc5dpp* ivna*(Ama1N + Ama1PN) ) * Cdc5AN) 
    dCdc5T_dt = 1.0*(kCdc5s*vna + kCdc5sp * Ndt80N - (kCdc5d + kCdc5dp * ivna*Ama1N + kCdc5dpp * ivna*(Ama1N + Ama1PN)) * Cdc5TN)
    dAma1P_dt = 1.0*((kAma1i + kAma1ip*ivna*ClbTN) * Ama1N*vna/(JAma1*vna + Ama1N) - kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - kAma1dp*Ama1PN)
    dAma1_dt = 1.0*(kAma1a * Ama1PN*vna/(JAma1*vna + Ama1PN) - (kAma1i + kAma1ip*ivna*ClbTN) * Ama1N*vna/(JAma1*vna + Ama1N)  + kAma1s*vna - kAma1dp*Ama1N + AlphakAma1s*(1*vna - Rim4(t,Rim4_tau, Rim4_T)*vna)*math.exp(-(t - Rim4_T)/25))
    
    # Return the result as a NumPy array
    return np.array([dClbT_dt, dNdt80_dt, dCdc20T_dt, dAPCp_dt, dAPCpCdc20_dt, dCdc5A_dt,  dCdc5T_dt, dAma1P_dt, dAma1_dt])

# Initial conditions
Odes_0 = np.array([(1.125* vna+0.12* vna+0.0),5 *vna,2*vna,0.0,0.1*vna,0.25 *vna,0.75*vna,0.0,0.0])

#Iteractive parameters definition 
pn.extension()
Rim4_tau_slider = pn.widgets.FloatSlider(
    name="Rim4 τ", start=0, end=10, step=0.01, value=0.2
)
Rim4_T_slider = pn.widgets.FloatSlider(
    name="Rim4 T", start=1, end=200, step=1, value=20
)
t_max_slider = pn.widgets.FloatSlider(
    name="t max",start=1, end=300, step=1, value=140
)
@pn.depends(
    Rim4_tau_slider.param.value,
    Rim4_T_slider.param.value,
    t_max_slider.param.value,
   
)

#Plotting function definition 
def ClbT_plot(
    Rim4_tau=0.2,
    Rim4_T=20,
    t_max=300,
):
    args = (
        Rim4_tau,
        Rim4_T
    )

      # Integrate ODES
    t = np.linspace(0, t_max, 400)
    Odes = scipy.integrate.odeint(ClbT_ODE, Odes_0, t, args=args)
    ClbTN, Ndt80N, Cdc20TN, APCpN, APCpCdc20N, Cdc5AN, Cdc5TN, Ama1PN, Ama1N  = Odes.transpose()
  
    # Set up plot
    p = bokeh.plotting.figure(
        frame_width=425,
        frame_height=350,
        x_axis_label="Time (minutes)",
        y_axis_label="Protein Concentration (a.u.)",
        x_range=[0, t_max],
        
    )

      # Plot of protein concentrations in arbitrary units 
    p.line(t, ClbTN/vna, line_width=2, color=colors[37], legend_label="ClbT")
    p.line(t, Ndt80N/vna, line_width=2, color=colors[3], legend_label="Ndt80N")
    p.line(t, Cdc20TN/vna, line_width=2, color=colors[4], legend_label="Cdc20TN")
    p.line(t, APCpN/vna, line_width=2, color=colors[87], legend_label="APCpN")
    p.line(t, APCpCdc20N/vna, line_width=2, color=colors[9], legend_label="APCpCdc20")
    p.line(t, Cdc5AN/vna, line_width=2, color=colors[70], legend_label="Cdc5A")
    p.line(t, Cdc5TN/vna, line_width=2, color=colors[8], legend_label="Cdc5TN")
    p.line(t, Ama1PN/vna, line_width=2, color=colors[14], legend_label="Ama1PN")
    p.line(t, Ama1N/vna , line_width=2, color=colors[10], legend_label="Ama1N")
    p.line(t,0.47,line_width=3,color=colors[11],line_dash='4 4')

    

    # Place the legend
    p.legend.location = "top_right"
    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"
    p.yaxis.major_label_text_font_size = '12pt'
    p.xaxis.major_label_text_font_size = '11pt'
    
    return p
    
widgets = pn.Column(
    pn.Spacer(height=10),
    Rim4_tau_slider,
    Rim4_T_slider,
    t_max_slider,
    width=150,
)

left_column = pn.Column(
    pn.Row(pn.Spacer(width=20), ClbT_plot,
))

# Final layout
pn.Row(left_column, pn.Spacer(width=20),widgets)


# In[ ]:




