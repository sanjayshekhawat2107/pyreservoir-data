"""
Material Balance Plots
@author: Sanjay Singh Shekhawat
@email: shekhawatsanjay2107@gmail.com
"""

def initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs):
  """
  Calculate OOIP and OGIP from Nfoi and Gfgi
  And output the result to labels in the plot
  """
  import matplotlib.patches as mpl_patches
  

  labels = []
  labels.append("Nfoi = {0:.4g} STB".format(Nfoi))
  labels.append("Gfgi = {0:.4g} SCF".format(Gfgi))
  labels.append("OOIP = {0:.4g} STB".format(OOIP))
  labels.append("OGIP = {0:.4g} SCF".format(OGIP))

  handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                  lw=0, alpha=0)] * 4
  return labels, handles, OOIP, OGIP   

class oil():
    """
    Oil (Undersaturated and saturated; Volatile and Non-volatile) Material Balance Plot
    """
    def calculate_params(self, p, Bo, Bg, Rv, Rs, Np, Gp, Gi, cf, cw, swi, Wp, We):
        """
        Calculate Material Balance Paramaters for Oil Reservoir
        
        Output: Fhc, F, Efw, Eo, Eg
        """
        pi = p[0]
        Rsi = Rs[0]
        Rvi = Rv[0]
        Boi = Bo[0]
        Bgi = Bg[0]

        # calculate Efw
        Efw = ((cf + cw * swi) / (1 - swi)) * (pi - p)

        # calculate Fhc, F
        Fhc = (Np * Bo) + (Gp-Gi) * (Rv-Rs)
        F = (Np * Bo) + (Gp-Gi) * (Rv-Rs) + Wp
        

        # calculate Eo and Eg
        Eo = Bo - Boi
        Eg = (Rsi-Rs)*Bg

        return Fhc, F, Efw, Eo, Eg

    def gascap(self, Gfgi, Nfoi, Bg, Bo):
      """
      Calculate Total Oil+Gas Expansion Factor from known Gas Cap ratio
      Gfgi and Nfoi known from volumetrics
      """
      Bgi, Boi = Bg[0], Bo[0]

      m = (Gfgi * Bgi) / (Nfoi * Boi)
      return m

    def plot(self, oil_type, F, Efw, Eo, Eg, Np, Bo, Rs, Rv, start=0, end=-1, figsize=(10,5)):
      """
      Create Material Balance Plots for Oil Reservoir
      
      Input:
      oil_type: 'saturated'
      """
      import numpy as np
      import matplotlib.pyplot as plt
      from scipy.optimize import curve_fit
      import matplotlib.patches as mpl_patches

      # plot attributes
      title_size = 15
      title_pad = 14

      # linear function for curve-fit
      def linear_zero_intercept(x, m):
          y = m * x
          return y

      def linear_with_intercept(x, m, c):
          y = m * x + c
          return y

      if oil_type == 'saturated':

        plt.figure(figsize=figsize)

        " Plot 1: Fhc/Eo vs Eg/Eo "

        plt.subplot(1,3,1)
        x1, y1 = (Eg / Eo), (Fhc / Eo)
        plt.plot(x1, y1, '.-')
        plt.title('Plot 1: Fhc/Eo vs Eg/Eo', size=title_size, pad=title_pad)
        plt.xlabel(r'$\frac{Eg}{Eo}$ (STB/scf)', size=15)
        plt.ylabel(r'$\frac{Fhc}{Eo}$ (STB)', size=15)

        ## curve-fitting to calculate the slope as Gfgi, intercept as Nfoi
        x1_norm = x1[1:] / max(x1[1:]) # normalize x
        y1_norm = y1[1:] / max(y1[1:]) # normalize y
        popt, pcov = curve_fit(linear_with_intercept, x1_norm, y1_norm)

        m, c = popt[0], popt[1]
        Gfgi = m = m * max(y1[1:]) / max(x1[1:]) # denormalize the slope
        Nfoi = c = c * max(y1[1:]) # denormalize the intercept

        ## calculate OOIP and OGIP from Nfoi and Gfgi
        Rsi, Rvi = Rs[0], Rv[0]
        OOIP = Nfoi + Gfgi * Rvi
        OGIP = Gfgi + Nfoi * Rsi

        ## Output results into text in plot
        labels, handles, OOIP, OGIP = initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs)    

        ## plot the regression line
        x1_fit = np.linspace(min(x1[1:]), max(x1[1:]), 5)
        y1_fit = linear_with_intercept(x1_fit, m, c)
        plt.plot(x1_fit, y1_fit)

        plt.legend(handles, labels, loc='best', fontsize='small', 
                   fancybox=True, framealpha=0.7, 
                   handlelength=0, handletextpad=0)

        " Plot 2: Fhc/Eg vs Eo/Eg "

        plt.subplot(1,3,2)
        x2, y2 =  (Eo / Eg), (Fhc / Eg)
        plt.plot(x2, y2, '.-')
        plt.title('Plot 2: Fhc/Eg vs Eo/Eg', size=title_size, pad=title_pad)
        plt.xlabel(r'$\frac{Eo}{Eg}$ (scf/STB)', size=15)
        plt.ylabel(r'$\frac{Fhc}{Eg}$ (scf)', size=15)

        ## curve-fitting to calculate the slope as Nfoi, intercept as Gfgi
        x2_norm = x2[1:] / max(x2[1:]) # normalize x
        y2_norm = y2[1:] / max(y2[1:]) # normalize y
        popt, pcov = curve_fit(linear_with_intercept, x2_norm, y2_norm)

        m, c = popt[0], popt[1]
        Nfoi = m = m * max(y2[1:]) / max(x2[1:]) # denormalize the slope
        Gfgi = c = c * max(y2[1:]) # denormalize the intercept

        ## calculate OOIP and OGIP from Nfoi and Gfgi
        Rsi, Rvi = Rs[0], Rv[0]
        OOIP = Nfoi + Gfgi * Rvi
        OGIP = Gfgi + Nfoi * Rsi

        ## Output results into text in plot
        labels, handles, OOIP, OGIP = initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs)    

        ## plot the regression line
        x2_fit = np.linspace(min(x2[1:]), max(x2[1:]), 5)
        y2_fit = linear_with_intercept(x2_fit, m, c)
        plt.plot(x2_fit, y2_fit)

        plt.legend(handles, labels, loc='best', fontsize='small', 
                   fancybox=True, framealpha=0.7, 
                   handlelength=0, handletextpad=0) 

        plt.tight_layout(1)                 

        plt.show()

        " Plot 3: F/(Eo + Eg) vs We/(Eo + Eg) "

        plt.subplot(1,3,2)
        x2, y2 =  F/(Eo+Eg), We/(Eo+Eg)
        plt.plot(x2, y2, '.-')
        plt.title('Plot 2: F/(Eo+Eg) vs We/(Eo+Eg)', size=title_size, pad=title_pad)
        plt.xlabel(r'$\frac{We}{Eo}{Eg}$ (STB)', size=15)
        plt.ylabel(r'$\frac{F}{Eo}{Eg}$ (STB)', size=15)

        ## curve-fitting to calculate the slope as Nfoi, intercept as Gfgi
        x2_norm = x2[1:] / max(x2[1:]) # normalize x
        y2_norm = y2[1:] / max(y2[1:]) # normalize y
        popt, pcov = curve_fit(linear_with_intercept, x2_norm, y2_norm)

        c = popt[0], popt[1]
        N = c = c * max(y2[1:]) # denormalize the intercept
  

        ## plot the regression line
        x2_fit = np.linspace(min(x2[1:]), max(x2[1:]), 5)
        y2_fit = linear_with_intercept(x2_fit, c)
        plt.plot(x2_fit, y2_fit)

        plt.legend(handles, labels, loc='best', fontsize='small', 
                   fancybox=True, framealpha=0.7, 
                   handlelength=0, handletextpad=0) 

        plt.tight_layout(1)                 

        plt.show()    

#     def plot(self, oil_type, F, Bto, Btg, Efw, Eo, Eg, Np, Bo, Rs, Rv, figsize=(10,5)):
#       """
#       Create Material Balance Plots for Oil Reservoir
      
#       Input:
#       oil_type: 'undersaturated' or 'saturated'
#       """
#       import numpy as np
#       import matplotlib.pyplot as plt
#       from scipy.optimize import curve_fit
#       import matplotlib.patches as mpl_patches

#       # plot attributes
#       title_size = 15
#       title_pad = 14

#       # linear function for curve-fit
#       def linear_zero_intercept(x, m):
#           y = m * x
#           return y

#       def linear_with_intercept(x, m, c):
#           y = m * x + c
#           return y

#       if oil_type == 'undersaturated':

#         plt.figure(figsize=figsize)

#         " Plot 1: F vs (Eg+Boi*Efw) "

#         plt.subplot(1,2,1)
#         Boi = Bo[0]
#         x1, y1 = (Eg + Boi * Efw), F
#         plt.plot(x1, y1, '.-')
#         plt.title(r'Plot 1: $F$ vs $(E_o+B_{oi}*E_{fw})$', size=title_size, pad=title_pad)
#         plt.xlabel(r'$E_o+B_{oi}E_{fw}$ (RB/STB)', size=15)
#         plt.ylabel(r'$F$ (res bbl)', size=15)

#         ## curve-fitting to calculate the slope as OOIP
#         x1_norm = x1 / max(x1) # normalize x
#         y1_norm = y1 / max(y1) # normalize y
#         popt, pcov = curve_fit(linear_zero_intercept, x1_norm, y1_norm)

#         m = popt[0]
#         Nfoi = m * max(y1) / max(x1) # denormalize the slope, hence the OGIP

#         ## Calculate OOIP and OGIP from Nfoi
#         Rsi = Rs[0]
#         Gfgi = 0 # no free gas phase in undersaturated oil
#         OOIP = Nfoi
#         OGIP = Nfoi * Rsi

#         ## Output results into text in plot
#         labels, handles, OOIP, OGIP = initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs)    

#         ## plot the regression line
#         x1_fit = np.linspace(min(x1), max(x1), 5)
#         y1_fit = linear_zero_intercept(x1_fit, Nfoi)
#         plt.plot(x1_fit, y1_fit, label='{} MMSTB'.format(np.round(Nfoi * 1E-6, 3)))

#         plt.legend(handles, labels, loc='best', fontsize='small', 
#                    fancybox=True, framealpha=0.7, 
#                    handlelength=0, handletextpad=0) 

#         " Plot 2: F/(Eg+Boi*Efw) vs Np (Waterdrive Diagnostic Plot) "

#         plt.subplot(1,2,2)
#         x2, y2 = Np, F / (Eg + Boi * Efw)
#         plt.plot(x2, y2, '.-')
#         plt.title('Plot 2: Waterdrive Diagnostic Plot', size=title_size, pad=title_pad)
#         plt.xlabel(r'$N_p$ (STB)', size=15)
#         plt.ylabel(r'$\frac{F}{(E_o+B_{oi}E_{fw})}$ (STB)', size=15)

#         ## curve-fitting to calculate the slope as OOIP, here [1:] because NaN is removed
#         x2_norm = x2[1:] / max(x2[1:]) # normalize x
#         y2_norm = y2[1:] / max(y2[1:]) # normalize y
#         popt, pcov = curve_fit(linear_with_intercept, x2_norm, y2_norm)

#         m, c = popt[0], popt[1]
#         m = m * max(y2[1:]) / max(x2[1:]) # denormalize the slope
#         Nfoi = c * max(y2[1:]) # denormalize the intercept, hence the OGIP

#         ## Calculate OOIP and OGIP from Nfoi
#         Rsi = Rs[0]
#         Gfgi = 0 # no free gas phase in undersaturated oil
#         OOIP = Nfoi
#         OGIP = Nfoi * Rsi

#         ## Output results into text in plot
#         labels, handles, OOIP, OGIP = initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs)           

#         ## plot the regression line
#         x2_fit = np.linspace(min(x2[1:]), max(x2[1:]), 5)
#         y2_fit = linear_with_intercept(x2_fit, m, Nfoi)
#         plt.plot(x2_fit, y2_fit, label='{} MMSTB'.format(np.round(Nfoi * 1E-6, 3)))

#         plt.legend(handles, labels, loc='best', fontsize='small', 
#                    fancybox=True, framealpha=0.7, 
#                    handlelength=0, handletextpad=0)  
        
#         plt.tight_layout(1)
#         plt.show()

#       if oil_type == 'saturated':

#         plt.figure(figsize=figsize)

#         " Plot 1: F/Eo vs Eg/Eo "

#         plt.subplot(1,3,1)
#         x1, y1 = (Eg / Eo), (F / Eo)
#         plt.plot(x1, y1, '.-')
#         plt.title('Plot 1: F/Eo vs Eg/Eo', size=title_size, pad=title_pad)
#         plt.xlabel(r'$\frac{Eg}{Eo}$ (STB/scf)', size=15)
#         plt.ylabel(r'$\frac{F}{Eo}$ (STB)', size=15)

#         ## curve-fitting to calculate the slope as Gfgi, intercept as Nfoi
#         x1_norm = x1[1:] / max(x1[1:]) # normalize x
#         y1_norm = y1[1:] / max(y1[1:]) # normalize y
#         popt, pcov = curve_fit(linear_with_intercept, x1_norm, y1_norm)

#         m, c = popt[0], popt[1]
#         Gfgi = m = m * max(y1[1:]) / max(x1[1:]) # denormalize the slope
#         Nfoi = c = c * max(y1[1:]) # denormalize the intercept

#         ## calculate OOIP and OGIP from Nfoi and Gfgi
#         Rsi, Rvi = Rs[0], Rv[0]
#         OOIP = Nfoi + Gfgi * Rvi
#         OGIP = Gfgi + Nfoi * Rsi

#         ## Output results into text in plot
#         labels, handles, OOIP, OGIP = initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs)    

#         ## plot the regression line
#         x1_fit = np.linspace(min(x1[1:]), max(x1[1:]), 5)
#         y1_fit = linear_with_intercept(x1_fit, m, c)
#         plt.plot(x1_fit, y1_fit)

#         plt.legend(handles, labels, loc='best', fontsize='small', 
#                    fancybox=True, framealpha=0.7, 
#                    handlelength=0, handletextpad=0)

#         " Plot 2: p/z vs Gp "

#         plt.subplot(1,3,2)
#         x2, y2 =  (Eo / Eg), (F / Eg)
#         plt.plot(x2, y2, '.-')
#         plt.title('Plot 2: F/Eg vs Eo/Eg', size=title_size, pad=title_pad)
#         plt.xlabel(r'$\frac{Eo}{Eg}$ (scf/STB)', size=15)
#         plt.ylabel(r'$\frac{F}{Eg}$ (scf)', size=15)

#         ## curve-fitting to calculate the slope as Nfoi, intercept as Gfgi
#         x2_norm = x2[1:] / max(x2[1:]) # normalize x
#         y2_norm = y2[1:] / max(y2[1:]) # normalize y
#         popt, pcov = curve_fit(linear_with_intercept, x2_norm, y2_norm)

#         m, c = popt[0], popt[1]
#         Nfoi = m = m * max(y2[1:]) / max(x2[1:]) # denormalize the slope
#         Gfgi = c = c * max(y2[1:]) # denormalize the intercept

#         ## calculate OOIP and OGIP from Nfoi and Gfgi
#         Rsi, Rvi = Rs[0], Rv[0]
#         OOIP = Nfoi + Gfgi * Rvi
#         OGIP = Gfgi + Nfoi * Rsi

#         ## Output results into text in plot
#         labels, handles, OOIP, OGIP = initial_hydrocarbon_in_place(Nfoi, Gfgi, Rv, Rs)    

#         ## plot the regression line
#         x2_fit = np.linspace(min(x2[1:]), max(x2[1:]), 5)
#         y2_fit = linear_with_intercept(x2_fit, m, c)
#         plt.plot(x2_fit, y2_fit)

#         plt.legend(handles, labels, loc='best', fontsize='small', 
#                    fancybox=True, framealpha=0.7, 
#                    handlelength=0, handletextpad=0) 

#         plt.tight_layout(1)                 

#         plt.show()    
