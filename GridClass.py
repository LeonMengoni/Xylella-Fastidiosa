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
                raise ValueError("Shape must be provided")
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

    def __set_control_zone(self): 
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
        
        self.density[self.EZ_mask] = 0 # eradicate all trees in EZ
        self.grove_mask = self.density > 0
        self.no_grove_mask = self.density == 0

    def __set_short_distance_kernel(self):
        y = np.linspace(-(self.rows-1), (self.rows-1), 2*self.rows-1)
        x = np.linspace(-(self.cols-1), (self.cols-1), 2*self.cols-1)
        Y, X = np.meshgrid(y, x)
        
        if self.kernel_type == "exponential":
            self.kernel = np.exp(-(X**2 + Y**2)**(1/2) / self.beta)

        elif self.kernel_type == "gaussian":
            self.kernel = np.exp(-(X**2 + Y**2)/(2 * self.beta**2)) / np.sqrt(2 * np.pi * self.beta**2)

    def __adjust_population(self, I):
        I[self.no_grove_mask] = 0
        I[I < self.tol] = 0
        return I

    # LOCAL GROWTH
    def __Gompertz_local_growth(self, I):
        I = self.K ** (1 - np.exp(-self.A)) * (I ** np.exp(-self.A)) # Number of locally infected trees
        I = self.__adjust_population(I)
        return I

    # SHORT DISTANCE DISPERSAL
    def __short_distance_dispersal(self, I):
        I = sp.signal.convolve(I, self.kernel, mode="same", method="fft")
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
                        I[tuple(new_cell)] += (1 - I[tuple(new_cell)]) * np.exp(-self.B) # Added the 1-I part to say that the remaining susceptible are infected with probability exp(-B) 
                        break    # OK
        I = self.__adjust_population(I)
        return I

    # LEVY FLIGHT DISPERSAL
    def __levy_flight_dispersal(self, I, V):
        if self.d_max is None:
            self.d_max = np.sqrt(np.floor(self.rows / 2)**2 + np.floor(self.cols / 2)**2); # about half the diagonal
        
        sample_power_law_params = {'alpha': self.alpha, 'x_min': self.d_min, 'x_max': self.d_max, 'sample': self.sample}
        d_list = ut.sample_power_law(size=V, **sample_power_law_params)

        for i, vector_pos, d in zip(range(self.V), self.vector_positions, d_list):
            while d > self.d_max: 
                d = ut.sample_power_law(size=1, **sample_power_law_params)

            while True:
                theta = np.random.uniform() * 2 * np.pi # random uniform direction
                step_coord = np.rint(d * np.array([np.sin(theta), np.cos(theta)])).astype(np.int32)
                new_vector_pos = vector_pos + step_coord
                if np.any(new_vector_pos < 0) or np.any(np.array(self.shape - new_vector_pos <= 0)): # Outside grid, REDO
                    continue
                elif self.sea_mask[tuple(new_vector_pos)]: # Sea, REDO
                    continue
                elif np.array_equal(new_vector_pos, vector_pos): # Same cell, OK
                    break
                elif self.no_grove_mask[tuple(new_vector_pos)]: # No grove, OK
                    break
                else:
                    I[tuple(new_vector_pos)] += (1 - I[tuple(new_vector_pos)]) * np.exp(-self.B) # Added the 1-I part to say that the remaining susceptible are infected with probability exp(-B) 
                    break
            self.vector_positions[i] = new_vector_pos
        I = self.__adjust_population(I)
        return I

    def __set_control_measures(self, I):
        rnd = np.random.random(size=self.shape)
        control_mask = (rnd < self.BZ_eff) & self.BZ_mask # & ~self.sea_mask if want to check over whole grid
        self.density[control_mask] = np.maximum(self.density[control_mask] - I[control_mask], 0) # Element-wise maximum of array elements
        I[control_mask] = 0
        self.grove_mask = self.density > 0 # update grove mask
        self.no_grove_mask = self.density == 0 # update no grove mask
        self.K = self.density + self.a * (1 - self.density) # update carrying capacity
        self.K[self.sea_mask] = 0
        return I
    
    def simulate(self, timesteps, parameters):
        self.timesteps = timesteps
        self.parameters = parameters

        ################################################################################
        # UNPACK PARAMETERS AND SET UP
        # Control Zone (CZ) parameters
        self.control, self.EZW, self.BZW, self.BZ_eff = self.parameters['control_zone']
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

        # Levy flight jumps
        elif self.dispersal_type == 'levy_flight':
            self.n_vectors, self.d_min, self.d_max, self.alpha, self.sample, self.vector_pop = parameters[self.dispersal_type]
        
            # Vector population dynamics (only if Levy flight dispersal)
            if self.vector_pop:
                self.max_vectors_per_cell, self.birth_rate, self.death_rate = parameters['vector_pop_dynamics']

        ################################################################################
        # DEFINE DATA ARRAYS
        self.I = np.zeros((self.timesteps+1, self.rows, self.cols)) # Number of infected trees (absolute density, maximum for a generic cell is self.density for that cell)
        self.incidence = np.zeros((self.timesteps+1, self.rows, self.cols)) # Fraction of infected trees
        
        if self.dispersal_type == 'levy_flight':
            self.V = np.zeros(self.timesteps+1) # vector population for all times
            
            if self.vector_pop:
                self.C = np.zeros(self.timesteps+1) # vector capacity (maximum number of vectors on the grid)

        ################################################################################
        # SET INITIAL (t=0) CONDITIONS
        self.I[0][tuple(self.seed)] = self.K[tuple(self.seed)] * np.exp(-self.B)
        self.incidence[0] = self.__evaluate_incidence(self.I[0])

        if self.dispersal_type == 'levy_flight':
            if not self.vector_pop:
                self.vector_positions = np.full((self.n_vectors,2), self.seed, dtype=int) # all vectors start in the seed cell
                self.V[0] = self.n_vectors
            if self.vector_pop:
                self.flag = False # used for initiating vector expansion (as soon as capacity is greater than 1)
                self.C[0] = self.__calculate_vector_capacity(self.I[0])
 
        ################################################################################
        # INITIATE EVOLUTION 
        # if self.vector_pop:
        #     self.flag = False

        for t in range(1, self.timesteps+1):
            # Local growth
            self.I[t] = self.__Gompertz_local_growth(self.I[t-1])

            # Dispersal
            if self.dispersal_type == 'short_long':
                self.I[t] = self.__short_distance_dispersal(self.I[t-1])
                self.I[t] = self.__long_distance_dispersal(self.I[t-1])

            elif self.dispersal_type == 'levy_flight':
                # print(self.V[t-1])
                self.I[t] = self.__levy_flight_dispersal(self.I[t-1], self.V[t-1])

                # Adjust vector population after the evolution
                if not self.vector_pop:
                    self.V[t] = self.n_vectors

                elif self.vector_pop:
                    self.C[t] = self.__calculate_vector_capacity(self.I[t])
                    self.V[t] = self.__adjust_vector_pop(self.V[t-1], self.C[t]) 
                    
                    # TODO: update vector positions, add or remove
                    # self.vector_positions = self.adjust_vector_positions()
                    new_vectors = self.V[t] - self.V[t-1]
                    if new_vectors > 0:
                        new_positions = self.vector_positions[np.random.choice(self.vector_positions.shape[0], size=new_vectors, replace=False)]
                        self.vector_positions = np.vstack((self.vector_positions, new_positions))
                    elif new_vectors == 0:
                        pass
                    else: 
                        remove_vectors = abs(new_vectors)
                        remove_indices = np.random.choice(self.vector_positions.shape[0], size=remove_vectors, replace=False)
                        mask = np.ones(self.vector_positions.shape[0], dtype=bool)
                        mask[remove_indices] = False
                        self.vector_positions = self.vector_positions[mask]

            # Control measures in buffer zone
            if self.from_file and self.control: 
                self.I[t] = self.__set_control_measures(self.I[t])
                
            # Obtain incidence
            self.incidence[t] = self.__evaluate_incidence(self.I[t])

    def __adjust_vector_pop(self, V, C):
        if self.flag:
            V = V + np.rint(self.birth_rate * V * (1 - V / self.C) - self.death_rate * self.V).astype(np.int32)
        if C > 1 and not self.flag:
            V = 1
            self.flag = True
        return V

    def __evaluate_incidence(self, I):
        incidence = np.zeros(I.shape)
        incidence[self.grove_mask] = I[self.grove_mask] / self.density[self.grove_mask]
        incidence[self.sea_mask] = -9999
        incidence[self.no_grove_mask] = 0
        return incidence
    
    def __calculate_vector_capacity(self, I):
        C = np.rint((I[self.grove_mask] * self.max_vectors_per_cell).sum()).astype(np.int32)
        return C

    def evaluate_risk(self, N, timesteps, parameters): # average incidence over N simulations
        self.N = N
        self.timesteps = timesteps
        self.parameters = parameters

        self.all_incidences = np.zeros((self.N, self.timesteps+1, self.rows, self.cols))
        # self.risk = np.zeros((self.timesteps+1, self.rows, self.cols))
        # self.final_incidence = np.zeros((self.N, self.rows, self.cols))

        # for i in range(N):
        #     self.simulate(self.timesteps, self.parameters)
        #     self.final_incidence[i] = self.incidence[-1]

        for i in range(N):
            self.simulate(self.timesteps, self.parameters)
            self.all_incidences[i] = self.incidence
        
        self.risk = np.mean(self.all_incidences, axis=0)

    def plot_short_distance_kernel(self, figsize=(6,6)):
        self.k_threshold = 1e-40
        self.k_min = max(self.k_threshold, np.min(self.kernel))
        logNorm = colors.LogNorm(vmin=self.k_min, vmax=1)

        fig, ax = plt.subplots(figsize=figsize)
        im = ax.imshow(self.kernel, norm=logNorm)
        if self.kernel_type == "exponential":   ax.set_title("Exponential short distance kernel")
        elif self.kernel_type == "gaussian":    ax.set_title("Gaussian short distance kernel")
        fig.colorbar(im, ax=ax)
        plt.show()

    def plot_density(self, figsize=(6,6)):
        cmap_custom = colors.LinearSegmentedColormap.from_list("", ["yellow", "forestgreen", "darkgreen"]) # grove density color map

        im = np.ma.array(self.density, mask=self.sea_mask)
        im_sea = np.ma.array(self.density, mask=~self.sea_mask)
        im_no_grove = np.ma.array(self.density, mask=~self.no_grove_mask)
        
        im_list = [im, im_sea, im_no_grove]
        color_list = ['tab:blue', 'black']
        cmap_list = [cmap_custom] + [colors.ListedColormap([color]) for color in color_list]
        label_list = ["Sea", "No groves"]

        if self.control and self.from_file:
            im_EZ = np.ma.array(self.density, mask=~self.EZ_mask)
            im_BZ = np.ma.array(self.density, mask=~self.BZ_mask)
            im_list += [im_EZ, im_BZ]
            new_colors = ['red', 'orange']
            color_list += new_colors
            cmap_list += [colors.ListedColormap([color]) for color in new_colors]
            label_list += ["Eradication zone (EZ)", "Buffer zone (BZ)"]

        fig, ax = plt.subplots(figsize=figsize)
        images = []
        for im, cmap in zip(im_list, cmap_list):
            images.append(ax.imshow(im, cmap=cmap, interpolation=None))

        ax.set_title("Olive Groves Density")
        cbar = fig.colorbar(images[0], ax=ax, orientation="horizontal")
        cbar.set_label("Grove density")
        patches = [mpatches.Patch(color=color_list[i], label=label_list[i]) for i in range(len(im_list)-1)]
        plt.legend(handles=patches, loc="upper right")

        plt.show()

    def plot_incidence(self, figsize=(6,6)):
        # Attributes for plot
        cmap_inferno = mpl.colormaps["inferno"]
        cmap_inferno.set_under("tab:blue")

        for t in range(len(self.incidence)):
            # Plot incidence for all times starting at 0
            fig, ax = plt.subplots(figsize=figsize)
            ax.set_title(f"Timestep {t}")
            im = ax.imshow(self.incidence[t], cmap=cmap_inferno, norm=colors.Normalize(vmin=0, vmax=1))
            cbar = fig.colorbar(im, ax=ax, orientation="horizontal")
            cbar.set_label("Disease incidence")
            plt.show()

    def plot_final_incidence(self, figsize=(6,6)):
        # Attributes for plot
        cmap_inferno = mpl.colormaps["inferno"]
        im_sea = np.ma.array(self.density, mask=~self.sea_mask)

        fig, ax = plt.subplots(figsize=figsize)
        ax.set_title(f"Final timestep {self.timesteps}")
        im = ax.imshow(self.incidence[-1], cmap=cmap_inferno, norm=colors.Normalize(vmin=0, vmax=1))
        ax.imshow(im_sea, cmap=colors.ListedColormap(['tab:blue']), interpolation=None)
        cbar = fig.colorbar(im, ax=ax, orientation="horizontal")
        cbar.set_label("Disease incidence")
        plt.show()