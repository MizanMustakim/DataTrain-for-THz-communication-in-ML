import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from scipy.integrate import quad


class IntraBody:
    def __init__(self, distance, frequency, ref_index, abs_coeff):
        self.distance = distance      # In mm
        self.frequency = frequency * 10**12             # In Hz
        self.ref_index = ref_index
        self.abs_coeff = abs_coeff
    
    def pathLoss(self):
        path_loss_dm = []
        for i in range(len(self.distance)):
            spreading = np.square(4*np.pi*self.frequency*self.ref_index*self.distance[i]*10**(-3) / sc.c)
            absn = np.exp(self.abs_coeff * self.distance[i])
            path_loss = spreading * absn
            path_loss_dm.append((10*np.log10(path_loss)))
        return path_loss_dm
    
    def noise_power(self):
        k = sc.value(u"Boltzmann constant")        # Boltzman constant in J/K
        h = sc.value(u"Planck constant")           # Planck's Constant in J/Hz
        temp = 310   # in Kelvin
        planck_func = ((2*sc.pi*h*self.frequency**3)/((sc.c)**2))/(np.exp((h*self.frequency)/(k*temp)))
        noise = []
        for i in range(len(self.distance)):
            func = lambda  freq:(-k*temp* (1-np.exp(-self.abs_coeff* self.distance[i]))* (sc.c/(np.sqrt(4*sc.pi)*freq))**2)
            val, err = quad(func,0, planck_func)
            noise.append(10*np.log10(val))
        return noise

def main():
    print("\n\n\n----------Start----------")
    
    distance = np.arange(0.1,(2.5+0.1),0.1)
    
    skin_path = []
    skin_noise = []
    
    blood_path = []
    blood_noise = []
    
    fat_path = []
    fat_noise = []
    
    for i in range(3):
        frequency = float(input("\n\nPlease Enter the frequency: "))
        ref_index = float(input("\n\nPlease Enter the Refractive Index\naccording to your given frequency: "))
        abs_coeff = float(input("\n\nPlease Enter the Absorption Coefficient\naccording to your given frequency: "))
        b = IntraBody(distance, frequency, ref_index,abs_coeff)
        path_loss = b.pathLoss()
        noise_power = b.noise_power()
        
        ### Medium Confirmation
        print("\n\n----------Please Choose the medium according to your Data----------")
        medium = input("\n\nS for Skin tissue\nB for Blood tissue\nF for fat tissue: ")
        
        if medium.upper() == "S":
            for j in range(len(path_loss)):
                skin_path.append(path_loss[j])
                skin_noise.append(noise_power[j])

        elif medium.upper() == "B":
            for k in range(len(path_loss)):
                blood_noise.append(noise_power[k])
                blood_path.append(path_loss[k])

        elif medium.upper() == "F":
            for l in range(len(path_loss)):
                fat_noise.append(noise_power[l])
                fat_path.append(path_loss[l])

        else:
            raise Exception("\n\nSorry! You entered a wrong input.\nPlease be sure and try again.")
        
    ### Plotting the graph for Path loss
    plt.plot(distance*10**(-3), blood_path, 'o-', label='blood tissue')
    plt.plot(distance*10**(-3), skin_path, '*-', label='skin tissue')
    plt.plot(distance*10**(-3), fat_path, '+-', label='fat tissue')
    plt.xlabel("Distance in m")
    plt.ylabel("Path Loss in dB")
    plt.title("Path loss in dB vs. Distance for different human tissues.")
    plt.legend()
    plt.show()
    

    ### Plotting the graph for Noise power
    plt.plot(distance*10**(-3), blood_noise, 'o-', label='blood tissue')
    plt.plot(distance*10**(-3), skin_noise, '*-', label='skin tissue')
    plt.plot(distance*10**(-3), fat_noise, '+-', label='fat tissue')
    plt.xlabel("Distance in m")
    plt.ylabel("Noise power in dB")
    plt.title("Noise power in dB vs. Distance for different human tissues.")
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()
        