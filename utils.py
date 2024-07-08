import numpy as np

def sample_power_law(size, alpha, x_min=0.1, x_max=None, sample="power"):
    u = np.random.random(size)

    if sample == "power":
        x = (1 - u)**(1 / (1 - alpha)) * x_min

    elif sample == "trunc":
        if x_max is None:
            raise Exception("x_max must be provided")
        x = (u * x_max**(1 - alpha) + (1 - u) * x_min**(1 - alpha))**(1 / (1 - alpha))

    return x