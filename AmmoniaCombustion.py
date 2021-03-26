import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

class CombustionModel():
    def __init__(self, gas, temperature, pressure, npoints):
        self.temperature = temperature
        self.pressure = pressure
        self.gas = gas
        self.npoints = npoints

    def calculate_and_plot(self):

        self.n_air_min = 0
        self.n_air_max = 8

        self.ppoints = np.size(self.pressure) # this is the size of pressure vector

        #create indexes
        self.io2 = self.gas.species_index('O2')
        self.ih2 = self.gas.species_index('H2')
        self.in2 = self.gas.species_index('N2')
        self.ino = self.gas.species_index('NO')
        self.ino2 = self.gas.species_index('NO2')
        self.in2o = self.gas.species_index('N2O')
        self.ih2o = self.gas.species_index('H2O')
        self.inh3 = self.gas.species_index('NH3')

        #gasifying agent: AIR
        self.x = np.zeros(self.gas.n_species) # this is the species list in gas

    #for plot no of mole fractions in the products
        self.z_o2 = np.zeros ((np.size(self.pressure),self.npoints))
        self.z_nh3 = np.zeros ((np.size(self.pressure),self.npoints))
        self.z_h2 = np.zeros ((np.size(self.pressure),self.npoints))
        self.z_n2 = np.zeros ((np.size(self.pressure),self.npoints))
        self.z_no = np.zeros ((np.size(self.pressure),self.npoints))
        self.z_no2 = np.zeros ((np.size(self.pressure),self.npoints))
        self.z_n2o = np.zeros ((np.size(self.pressure),self.npoints))
        self.z_h2o = np.zeros ((np.size(self.pressure),self.npoints))

        self.t_adiabatic = np.zeros ((np.size(self.pressure),self.npoints))

        self.n_ratio = np.zeros(self.npoints)

        self.d_n =  (self.n_air_max - self.n_air_min)/self.npoints

        for j in range(0,self.ppoints):

            for i in range (0,self.npoints):

                self.n_o2 = 0.21*(self.n_air_min + i*self.d_n)
                self.n_nh3 = 1 #mol
                self.n_n2 =0.79*(self.n_air_min+i*self.d_n)
                self.n_air = self.n_o2+self.n_n2


                self.x[self.inh3] = self.n_nh3/(self.n_nh3+self.n_o2+self.n_n2)
                self.x[self.io2] = self.n_o2/(self.n_nh3+self.n_o2+self.n_n2)
                self.x[self.in2] = self.n_n2/(self.n_nh3+self.n_o2+self.n_n2)

                self.n_ratio[i] = self.n_air/self.n_nh3

                self.gas.TPX = self.temperature, self.pressure[j], self.x

                #equilibration
                self.gas.equilibrate('HP', solver='auto')

                # defining i position inside zero matrix
                self.z_o2[j][i]= self.gas.X[self.io2]
                self.z_h2[j][i]= self.gas.X[self.ih2]
                self.z_n2[j][i]= self.gas.X[self.in2]
                self.z_no[j][i]= self.gas.X[self.ino]
                self.z_no2[j][i] = self.gas.X[self.ino2]
                self.z_n2o[j][i] = self.gas.X[self.in2o]
                self.z_h2o[j][i] = self.gas.X[self.ih2o]
                self.z_nh3[j][i] = self.gas.X[self.inh3]

                # i position of variables
                self.t_adiabatic[j][i] = self.gas.T


                self.m = ['x','o','s','>']

        plt.figure(figsize=(8, 6))
        for j in range(0,self.ppoints):
            plt.plot((self.n_ratio), self.t_adiabatic[j][:], label = '%s MPa' %(self.pressure[j]/10**6), marker = self.m[j],markevery = 10)

        plt.grid()
        plt.xlabel('$n_{air}/n_{NH3}$, [mole/mole]')
        plt.ylabel('$T$, [K]')
        plt.legend(loc='center right')
        plt.savefig('T_pressures_air.png',dpi=300)


        fig = plt.figure(figsize=(8, 6),facecolor='w', edgecolor='k')
        for j in range(0,self.ppoints):
            ax = fig.add_subplot(2, 2,j+1)
            ax.plot(self.n_ratio, self.z_h2[j][:], label='H2', marker = ".", markevery = 20)
            ax.plot(self.n_ratio, self.z_o2[j][:], label='O2',marker = ">",markevery = 20)
            ax.plot(self.n_ratio, self.z_h2o[j][:], label='H2O',marker = "+",markevery = 20)
            ax.plot(self.n_ratio, self.z_n2[j][:], label='N2',marker = "*",markevery = 20)
            ax.plot(self.n_ratio, self.z_nh3[j][:], label='NH3', marker ="x",markevery = 20)

            ax.grid()
            ax.legend(fancybox=True, framealpha=0.5, loc = 'upper right')
            if j>1:
                ax.set_xlabel('$n_{air}/n_{NH3}$, [mole/mole]')

            if (j==0 or j==2):
                ax.set_ylabel('Mole fraction $x_i$, [mole/mole]')

            tittext = '$p$ = ' + str(self.pressure[j]/10**6) + ' [MPa]'
            plt.title(tittext,fontsize = 10)

        #plt.show()
        plt.savefig('species_presssures_air.png',dpi=300)


        fig_1 = plt.figure(figsize=(8, 6),facecolor='w', edgecolor='k')
        for j in range(0,self.ppoints):
            ax = fig_1.add_subplot(2, 2,j+1)
            ax.set_yscale('log')
            ax.set_ylim(10**(-20),0.1)
            ax.plot(self.n_ratio, self.z_no[j][:], label='NO',marker = "o",markevery = 18)
            ax.plot(self.n_ratio, self.z_no2[j][:], label='NO2',marker = "s",markevery = 20)
            ax.plot(self.n_ratio, self.z_n2o[j][:], label='N2O',marker = "x",markevery = 22)

            ax.grid()
            ax.legend(fancybox=True, framealpha=0.5, loc = 'lower right')
            if j>1:
                ax.set_xlabel('$n_{air}/n_{NH3}$, [mole/mole]')

            if (j==0 or j==2):
                ax.set_ylabel('Mole fraction $x_i$, [mole/mole]')

            tittext = '$p$ = ' + str(self.pressure[j]/10**6) + ' [MPa]'
            plt.title(tittext,fontsize = 10)

        #plt.show()
        plt.savefig('NOx_presssures_air.png',dpi=300)
        return

if __name__ == '__main__':
    Ammonia = CombustionModel(ct.Solution('chem_ammonia.cti'), 300, [10**5,10**6,5*10**6,10**7], 300)
    Ammonia.calculate_and_plot()

