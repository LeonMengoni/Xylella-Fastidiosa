import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
import scipy as sp
import matplotlib.patches as mpatches
import utils as ut

class Grid():
    def __init__(self, density=None, shape=None, sea_perc=0.25):
        if density is not None: # from file
            self.from_file = True # set flag
            self.density = density
            self.shape = self.density.shape
            self.sea_mask = self.density < 0
            self.no_grove_mask = self.density == 0
            
        else: # generate
            if shape is None:
                raise ValueError('Shape must be provided')
            self.from_file = False # set flag
            self.shape = shape
            self.density = np.random.random(self.shape)

            # set random sea cells (density < 0)
            self.sea_mask = np.random.choice([True, False], size=self.shape, p=[sea_perc, 1-sea_perc])
            self.density[self.sea_mask] = -9999

            # set random null cells (olive_prop = 0)
            zero_perc = 0.1 # including sea
            zero_mask = np.random.choice([True, False], size=self.shape, p=[zero_perc, 1-zero_perc])
            self.no_grove_mask = np.logical_and(~self.sea_mask, zero_mask)
            self.density[self.no_grove_mask] = 0

        self.rows, self.cols = self.shape
        self.grove_mask = self.density > 0
        self.control = False # control zone flag

        # Set initial infection seed
        if self.from_file:
            self.seed = np.array([235,266]) # Initial spread from Gallipoli
        else:
            grove_coordinates = np.argwhere(self.grove_mask)
            self.seed = grove_coordinates[np.random.choice(grove_coordinates.shape[0])] # Random seed

    def __set_control_zone(self): # define the EZ and BZ
        xp1, yp1 = [253, 212]
        xp2, yp2 = [278, 183]

        m = (yp1 - yp2) / (xp1 - xp2)
        B1 = yp1 - m * xp1 # y-intercept of front-line/start of control zone
        B2_erad = B1 - self.EZW * np.sqrt(m**2 + 1) # y-intercept of eradication zone limit
        B2_rog = B1 - (self.EZW + self.BZW) * np.sqrt(m**2 + 1) # y-intercept of buffer zone limit/end of control zone

        self.EZ_mask = np.zeros(shape=self.shape, dtype=bool)
        self.BZ_mask = np.zeros(shape=self.shape, dtype=bool)
        for y in range(self.rows):
            for x in range(self.cols):
                if y < m * x + B1 and y > m * x + B2_erad and self.density[y,x] >= 0:
                    self.EZ_mask[y,x] = True
                elif y < m*x+B2_erad and y > m*x+B2_rog and self.density[y,x] >= 0:
                    self.BZ_mask[y,x] = True
        
        # self.density[self.EZ_mask] = 0 # eradicate all trees in EZ
        # self.grove_mask = self.density > 0
        # self.no_grove_mask = self.density == 0

    def __set_short_distance_kernel(self): # not normalized so that the value in (x,y)=(0,0) is 1
        y = np.linspace(-(self.rows-1), (self.rows-1), 2*self.rows-1)
        x = np.linspace(-(self.cols-1), (self.cols-1), 2*self.cols-1)
        Y, X = np.meshgrid(y, x)
        
        if self.kernel_type == 'exponential':
            self.kernel = np.exp(-(X**2 + Y**2)**(1/2) / self.beta)

        elif self.kernel_type == 'gaussian':
            self.kernel = np.exp(-(X**2 + Y**2)/(2 * self.beta**2))
        
    def __set_levy_flight_kernel(self):
        y = np.linspace(-(self.rows-1), (self.rows-1), 2*self.rows-1)
        x = np.linspace(-(self.cols-1), (self.cols-1), 2*self.cols-1)
        Y, X = np.meshgrid(y, x)
        Y[Y == 0] = 1
        X[X == 0] = 1

        # self.kernel = (1 - self.alpha) / (self.d_max**(1 - self.alpha) - self.d_min**(1 - self.alpha)) * (X**2 + Y**2)#**(-self.alpha / 2)
        self.kernel = (X**2 + Y**2)**(-self.alpha / 2)

    def __adjust_population(self, I):
        I[self.no_grove_mask] = 0
        I[I < self.tol] = 0
        return I

    # LOCAL GROWTH
    def __Gompertz_local_growth(self, I):
        I = self.K ** (1 - np.exp(-self.A)) * (I ** np.exp(-self.A)) # Number of locally infected trees
        I = self.__adjust_population(I)
        return I

    # KERNEL DISPERSAL
    def __kernel_convolution(self, I):
        I = sp.signal.convolve(I, self.kernel, mode='same', method='fft')
        I = self.__adjust_population(I)
        return I

    # LONG DISTANCE DISPERSAL
    def __long_distance_dispersal(self, I):
        prob_disp = np.random.random(size=self.shape) * I # dispersal probability for every cell
        disp_cells = np.argwhere(prob_disp > self.disp_tol) # cells that disperse
        rnd_disp = np.random.randint(1, self.M_max+1, size=len(disp_cells)) # random number of dispersers per cell

        for disp_cell, n_disp in zip(disp_cells, rnd_disp):
            for i_disp in range(n_disp):
                while True:
                    new_cell = disp_cell + np.rint(np.random.normal(0, self.D, size=2)).astype(np.int32)
                    if np.any(new_cell < 0) or np.any(np.array(self.shape - new_cell <= 0)):    continue # Outside grid
                    elif np.array_equal(new_cell, disp_cell):                                   continue # Same cell
                    elif self.sea_mask[tuple(new_cell)]:                                        continue # Sea
                    elif self.no_grove_mask[tuple(new_cell)]:                                   break    # No grove, OK
                    else:
                        I[tuple(new_cell)] += (self.density[tuple(new_cell)] - I[tuple(new_cell)]) * np.exp(-self.B) # Added the 1-I part to say that the remaining susceptible are infected with probability exp(-B) 
                        break    # OK
        I = self.__adjust_population(I)
        return I

    def __set_control_measures(self, I):
        rnd = np.random.random(size=self.shape)
        EZ_control_mask = (rnd < self.EZ_eff) & self.EZ_mask # if EZ_eff = 1, EZ_control_mask = EZ_mask
        BZ_control_mask = (rnd < self.BZ_eff) & self.BZ_mask # & ~self.sea_mask if want to check over whole grid
    
        self.density[EZ_control_mask] = np.maximum(self.density[EZ_control_mask] - I[EZ_control_mask], 0) # Element-wise maximum of array elements    
        self.density[BZ_control_mask] = np.maximum(self.density[BZ_control_mask] - I[BZ_control_mask], 0) 
        
        I[EZ_control_mask] = 0
        I[BZ_control_mask] = 0

        self.grove_mask = self.density > 0 # update grove mask
        self.no_grove_mask = self.density == 0 # update no grove mask
        self.K = self.density + self.a * (1 - self.density) # update carrying capacity
        self.K[self.sea_mask] = 0
        return I

    def __evaluate_incidence(self, I):
        incidence = np.zeros(I.shape)
        incidence[self.grove_mask] = I[self.grove_mask] / self.density[self.grove_mask]
        incidence[self.sea_mask] = -9999
        incidence[self.no_grove_mask] = 0
        return incidence

    # SIMULATE DIFFUSION
    def simulate(self, timesteps, parameters):
        self.timesteps = timesteps
        self.parameters = parameters
        
        ################################################################################
        # UNPACK PARAMETERS AND SET UP
        # Control Zone (CZ) parameters
        self.control, self.EZW, self.BZW, self.EZ_eff, self.BZ_eff = self.parameters['control_zone']
        if self.from_file and self.control:
            self.__set_control_zone()
        
        # Common parameters
        self.A, self.B, self.a, self.tol = self.parameters['common']
        self.K = self.density + self.a * (1 - self.density) # carrying capacity
        self.K[self.sea_mask] = 0

        # Dispersal type
        self.dispersal_type, = self.parameters['dispersal']

        # Short + long kernel
        if self.dispersal_type == 'short_long':
            self.beta, self.kernel_type, self.disp_tol, self.M_max, self.D = parameters[self.dispersal_type]
            self.__set_short_distance_kernel()

        elif self.dispersal_type == 'levy_flight':
            self.d_min, self.d_max, self.alpha = parameters[self.dispersal_type]
            self.__set_levy_flight_kernel()

        ################################################################################
        # DEFINE DATA ARRAYS
        self.I = np.zeros((self.timesteps+1, self.rows, self.cols)) # Number of infected trees (absolute density, maximum for a generic cell is self.density for that cell)
        self.incidence = self.I.copy() # Fraction of infected trees

        ################################################################################
        # SET INITIAL (t=0) CONDITIONS
        self.I[0][tuple(self.seed)] = self.K[tuple(self.seed)] * np.exp(-self.B)
        self.incidence[0] = self.__evaluate_incidence(self.I[0])

        ################################################################################
        # INITIATE EVOLUTION
        for t in range(1, self.timesteps + 1):
            # Local growth
            self.I[t] = self.__Gompertz_local_growth(self.I[t-1])

            if self.dispersal_type == 'short_long':
                self.I[t] = self.__kernel_convolution(self.I[t])
                self.I[t] = self.__long_distance_dispersal(self.I[t])
            
            elif self.dispersal_type == 'levy_flight':
                self.I[t] = self.__kernel_convolution(self.I[t])

            # Control efficiency in buffer zone
            if self.from_file and self.control: 
                self.I[t] = self.__set_control_measures(self.I[t])

            # Obtain incidence
            self.incidence[t] = self.__evaluate_incidence(self.I[t])

    def evaluate_risk(self, N, timesteps, parameters, verbose=True): # average incidence over N simulations
        self.N = N
        self.timesteps = timesteps
        self.parameters = parameters

        self.all_incidences = np.zeros((self.N, self.timesteps+1, self.rows, self.cols))

        for i in range(N):
            if verbose:
                print(f'Simulating run {i+1} of {N}...')
            self.simulate(self.timesteps, self.parameters)
            self.all_incidences[i] = self.incidence
        
        self.risk = np.mean(self.all_incidences, axis=0)

    def plot_short_distance_kernel(self, figsize=(6,6)):
        self.k_threshold = 1e-40
        self.k_min = max(self.k_threshold, np.min(self.kernel))
        logNorm = colors.LogNorm(vmin=self.k_min, vmax=1)

        fig, ax = plt.subplots(figsize=figsize)
        im = ax.imshow(self.kernel, norm=logNorm)
        if self.kernel_type == 'exponential':   ax.set_title('Exponential short distance kernel')
        elif self.kernel_type == 'gaussian':    ax.set_title('Gaussian short distance kernel')
        fig.colorbar(im, ax=ax)
        plt.show()

    def plot_density(self, figsize=(6,6)):
        cmap_custom = colors.LinearSegmentedColormap.from_list('', ['yellow', 'forestgreen', 'darkgreen']) # grove density color map

        im = np.ma.array(self.density, mask=self.sea_mask)
        im_sea = np.ma.array(self.density, mask=~self.sea_mask)
        im_no_grove = np.ma.array(self.density, mask=~self.no_grove_mask)
        
        im_list = [im, im_sea, im_no_grove]
        color_list = ['tab:blue', 'black']
        cmap_list = [cmap_custom] + [colors.ListedColormap([color]) for color in color_list]
        label_list = ['Sea', 'No groves']

        if self.control and self.from_file:
            im_EZ = np.ma.array(self.density, mask=~self.EZ_mask)
            im_BZ = np.ma.array(self.density, mask=~self.BZ_mask)
            im_list += [im_EZ, im_BZ]
            new_colors = ['red', 'orange']
            color_list += new_colors
            cmap_list += [colors.ListedColormap([color]) for color in new_colors]
            label_list += ['Eradication zone (EZ)', 'Buffer zone (BZ)']

        fig, ax = plt.subplots(figsize=figsize)
        images = []
        for im, cmap in zip(im_list, cmap_list):
            images.append(ax.imshow(im, cmap=cmap, interpolation=None))

        ax.set_title('Olive Groves Density')
        cbar = fig.colorbar(images[0], ax=ax, orientation='horizontal')
        cbar.set_label('Grove density')
        patches = [mpatches.Patch(color=color_list[i], label=label_list[i]) for i in range(len(im_list)-1)]
        plt.legend(handles=patches, loc='upper right')

        plt.show()

    def plot_incidence(self, zoom=False, figsize=(6,6)):
        # Attributes for plot
        cmap_inferno = mpl.colormaps['inferno']
        cmap_inferno.set_under('tab:blue')

        for t in range(len(self.incidence)):
            # Plot incidence for all times starting at 0
            fig, ax = plt.subplots(figsize=figsize)
            ax.set_title(f'Timestep {t}')
            if zoom:
                im = ax.imshow(self.incidence[t,150:,200:], cmap=cmap_inferno, norm=colors.Normalize(vmin=0, vmax=1))
            else: 
                im = ax.imshow(self.incidence[t], cmap=cmap_inferno, norm=colors.Normalize(vmin=0, vmax=1))
            cbar = fig.colorbar(im, ax=ax, orientation='horizontal')
            cbar.set_label('Disease incidence')
            plt.show()

    def plot_final_incidence(self, figsize=(6,6)):
        # Attributes for plot
        cmap_inferno = mpl.colormaps['inferno']
        im_sea = np.ma.array(self.density, mask=~self.sea_mask)

        fig, ax = plt.subplots(figsize=figsize)
        ax.set_title(f'Final timestep {self.timesteps}')
        im = ax.imshow(self.incidence[-1], cmap=cmap_inferno, norm=colors.Normalize(vmin=0, vmax=1))
        ax.imshow(im_sea, cmap=colors.ListedColormap(['tab:blue']), interpolation=None)
        cbar = fig.colorbar(im, ax=ax, orientation='horizontal')
        cbar.set_label('Disease incidence')
        plt.show()