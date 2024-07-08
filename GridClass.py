import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
import scipy as sp
import matplotlib.patches as mpatches
from utils import *

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

    def plot_density(self, figsize=(6,6)):
        cmap_custom = colors.LinearSegmentedColormap.from_list("", ["yellow", "forestgreen", "darkgreen"]) # grove density color map

        im = np.ma.array(self.density, mask=self.sea_mask)
        im_sea = np.ma.array(self.density, mask=~self.sea_mask)
        im_no_grove = np.ma.array(self.density, mask=~self.no_grove_mask)
        
        im_list = [im, im_sea, im_no_grove]
        color_list = ['tab:blue', 'black']
        cmap_list = [cmap_custom] + [colors.ListedColormap([color]) for color in color_list]
        label_list = ["Sea", "No groves"]

        if self.control:
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

    def set_control_zone(self, EZW, BZW, BZ_eff=0.9): 
        # EZW: eradication zone width
        # BZW: buffer zone width 
        # BZ_eff: buffer zone roguing efficiency
        self.control = True
        self.EZW = EZW
        self.BZW = BZW
        self.BZ_eff = BZ_eff
        xp1, yp1 = [253, 212]
        xp2, yp2 = [278, 183]

        m = (yp1 - yp2) / (xp1 - xp2)
        B1 = yp1 - m * xp1 # y-intercept of front-line/start of control zone
        B2_erad = B1 - EZW * np.sqrt(m**2 + 1) # y-intercept of eradication zone limit
        B2_rog = B1 - (EZW + BZW) * np.sqrt(m**2 + 1) # y-intercept of buffer zone limit/end of control zone

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

    def set_short_distance_kernel(self, beta, type="exponential"):
        self.beta = beta
        self.kernel_type = type

        y = np.linspace(-(self.rows-1), (self.rows-1), 2*self.rows-1)
        x = np.linspace(-(self.cols-1), (self.cols-1), 2*self.cols-1)
        Y, X = np.meshgrid(y, x)
        
        if type == "exponential":
            self.kernel = np.exp(-(X**2 + Y**2)**(1/2) / self.beta)

        elif type == "gaussian":
            self.kernel = np.exp(-(X**2 + Y**2)/(2 * self.beta**2)) / np.sqrt(2 * np.pi * self.beta**2)

    def plot_short_distance_kernel(self, figsize=(6,6)):
        k_threshold = 1e-40
        k_min = max(k_threshold, np.min(self.kernel))
        logNorm = colors.LogNorm(vmin=k_min, vmax=1)

        fig, ax = plt.subplots(figsize=figsize)
        im = ax.imshow(self.kernel, norm=logNorm)
        if self.kernel_type == "exponential":   ax.set_title("Exponential short distance kernel")
        elif self.kernel_type == "gaussian":    ax.set_title("Gaussian short distance kernel")
        fig.colorbar(im, ax=ax)
        plt.show()

    # LOCAL GROWTH
    def __Gompertz_local_growth(self):
        self.I = self.K ** (1 - np.exp(-self.A)) * (self.I ** np.exp(-self.A)) # Number of locally infected trees

    def __adjust_population(self):
        self.I[self.no_grove_mask] = 0
        self.I[self.I < self.tol] = 0

    # SHORT DISTANCE DISPERSAL
    def __short_distance_dispersal(self):
        self.I = sp.signal.convolve(self.I, self.kernel, mode="same", method="fft")

    # LONG DISTANCE DISPERSAL
    def __long_distance_dispersal(self):
        # self.long_distance_parameters = long_distance_parameters
        # disp_tol, M_max, D, B = self.long_distance_parameters

        prob_disp = np.random.random(size=self.shape) * self.I # dispersal probability for every cell
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
                        self.I[tuple(new_cell)] += (1 - self.I[tuple(new_cell)]) * np.exp(-self.B) # Added the 1-I part to say that the remaining susceptible are infected with probability exp(-B) 
                        break    # OK

    # LEVY FLIGHT DISPERSAL
    # TODO: add "scale" parameter for rescaling of random samples
    # TODO: Vectors have a certain probability of infecting a susceptible tree.
    # TODO: Vectors can be eliminated, according to a certain probability, which depends on preventive measures implemented: weeding, roguing, pesticides, etc.
    def __levy_flight_dispersal(self):
        if self.sample == "power":
            self.d_max = None
        elif self.sample == "trunc":
            self.d_max = np.sqrt(np.floor(self.rows / 2)**2 + np.floor(self.cols / 2)**2); # about half the diagonal
        
        sample_power_law_params = {'alpha': self.alpha, 'x_min': self.d_min, 'x_max': self.d_max, 'sample': self.sample}
        d_list = sample_power_law(size=self.n_vectors, **sample_power_law_params)

        for i, vector_pos, d in zip(range(self.n_vectors), self.vector_positions, d_list):
            y, x = vector_pos

            # If jump is larger than grid, generate new jump
            max_square_dist = np.max([x**2 + y**2,
                                    x**2 + (self.rows - y)**2,
                                    (self.cols - x)**2 + (self.rows - y)**2,
                                    (self.cols - x)**2 + y**2]) - 1 # -1 to account for discreteness of grid
            while d**2 > max_square_dist:
                d = sample_power_law(size=1, **sample_power_law_params)

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
                    self.I[tuple(new_vector_pos)] += (1 - self.I[tuple(new_vector_pos)]) * np.exp(-self.B) # Added the 1-I part to say that the remaining susceptible are infected with probability exp(-B) 
                    break
            self.vector_positions[i] = new_vector_pos

    # SIMULATE DIFFUSION
    def simulate(self, timesteps, parameters, dispersal_type="short_long", plot=True):
        # # Gompertz growth function parameters
        # A: rate of population growth
        # B: related to initial infections

        # a: relative carrying capacity in non olive plants
        # tol: numerical noise tolerance (IF TOLERANCE IS TOO BIG, IT CAN HAPPEN THAT THE SIMULATION IS TRUNCATED AT THE FIRST ITERATION)

        # # Short distance kernel
        # beta: mean dispersal distance

        # # Long distance kernel
        # disp_tol: stratified dispersal tolerance (threshold probability)
        # M_max: maximum number of dispersers
        # D: dispersal standard deviation

        # # Levy flight
        # n_vectors: number of disease-spreading vectors
        # d_min: minimum Levy jump
        # alpha: power law distribution exponent
        # sample: kwarg for type of power law distribution to sample from

        self.timesteps = timesteps
        self.parameters = parameters
        self.dispersal_type = dispersal_type
        
        self.output = np.zeros((self.timesteps+1, self.rows, self.cols)) # Fraction of infected trees
        self.I = np.zeros(self.shape) # Number of infected trees

        # Unpack common parameters
        self.A, self.B, self.a, self.tol = self.parameters['common']
        self.K = self.density + self.a * (1 - self.density) # carrying capacity
        self.K[self.sea_mask] = 0

        # Unpack parameters for different dispersal types
        if self.dispersal_type == 'short_long':
            self.beta, self.disp_tol, self.M_max, self.D = parameters[self.dispersal_type]
        elif self.dispersal_type == 'levy_flight':
            self.n_vectors, self.d_min, self.alpha, self.sample = parameters[self.dispersal_type]
            self.vector_positions = np.full((self.n_vectors,2), self.seed, dtype=int)
        
        self.I[tuple(self.seed)] = self.K[tuple(self.seed)] * np.exp(-self.B)

        # initiate evolution
        for t in range(self.output.shape[0]):
            if t > 0:
                # Local growth
                self.__Gompertz_local_growth()
                self.__adjust_population()

                # Levy flight dispersal
                if self.dispersal_type == 'levy_flight':
                    self.__levy_flight_dispersal()
                    self.__adjust_population()

                if self.dispersal_type == 'short_long':
                    # Short distance dispersal
                    self.__short_distance_dispersal()
                    self.__adjust_population()
                    
                    # Long distance dispersal
                    self.__long_distance_dispersal()
                    self.__adjust_population()
            
            # Control efficiency in buffer zone
            rnd = np.random.random(size=self.shape)
            control_mask = (rnd < self.BZ_eff) & self.BZ_mask # & ~self.sea_mask if want to check over whole region
            self.density[control_mask] = np.maximum(self.density[control_mask] - self.I[control_mask], 0) # Element-wise maximum of array elements
            self.I[control_mask] = 0
            self.grove_mask = self.density > 0 # update grove mask
            self.no_grove_mask = self.density == 0 # update no grove mask
            self.K = self.density + self.a * (1 - self.density) # update carrying capacity
            self.K[self.sea_mask] = 0

            # Obtain output
            self.output[t][self.grove_mask] = self.I[self.grove_mask] / self.density[self.grove_mask]
            self.output[t][self.sea_mask] = -9999
            self.output[t][self.no_grove_mask] = 0

            # if plot:
            #     self.plot_incidence(t)

    def plot_incidence(self, figsize=(6,6)):
        # Attributes for plot
        cmap_inferno = mpl.colormaps["inferno"]
        cmap_inferno.set_under("tab:blue")

        for t in range(len(self.output)):
            # Plot output for all times starting at 0
            fig, ax = plt.subplots(figsize=figsize)
            ax.set_title(f"Timestep {t}")
            im = ax.imshow(self.output[t], cmap=cmap_inferno, norm=colors.Normalize(vmin=0, vmax=1))
            cbar = fig.colorbar(im, ax=ax, orientation="horizontal")
            cbar.set_label("Disease incidence")
            plt.show()

    def plot_final_incidence(self, figsize=(6,6)):
        # Attributes for plot
        cmap_inferno = mpl.colormaps["inferno"]
        im_sea = np.ma.array(self.density, mask=~self.sea_mask)

        fig, ax = plt.subplots(figsize=figsize)
        ax.set_title(f"Final timestep {self.timesteps}")
        im = ax.imshow(self.output[-1], cmap=cmap_inferno, norm=colors.Normalize(vmin=0, vmax=1))
        ax.imshow(im_sea, cmap=colors.ListedColormap(['tab:blue']), interpolation=None)
        cbar = fig.colorbar(im, ax=ax, orientation="horizontal")
        cbar.set_label("Disease incidence")
        plt.show()