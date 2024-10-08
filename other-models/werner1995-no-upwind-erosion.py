import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys 
sys.setrecursionlimit(10_000)


'''
Werner 
+ with 100% deposition in shadow

- no upwind erosion, see fig 5b in https://smallpond.ca/jim/sand/dunefieldMorphology/index.html
'''



# Set up the simulation parameters
H = 5

width, height = 100, 100 
#grid = H*np.ones((height, width))  
grid = np.random.randint(3, 8, size=(width, height))
delta_avalanche = 5

wind_x, wind_y =  1, 0
p_sand = 0.6
p_surface = 0.4


shadow_coeff = np.tan(np.deg2rad(36))
print(shadow_coeff)

###########################################################################################################

def check_avalanche_condition(i, j):

    # no upwind erosion (upwind ~ not shadow region)
    if not is_in_shadow(i, j):
        return None
    
    l = (i+1)%width
    r = (i-1)%width
    u = (j+1)%height
    d = (j-1)%height
    bordering_cell_indizes = (l,j), (r,j), (i,u), (i,d)

    delta_heights = [grid[i,j] - grid[bc] for bc in bordering_cell_indizes]
    delta_max_index = np.argmax(delta_heights) # maximum height

    # call function to move sand particle from highest to lowest cell (that exceeds the delta)
    if delta_heights[delta_max_index] >= delta_avalanche:
        remove_particle(i, j) 
        add_particle(*bordering_cell_indizes[delta_max_index])

###########################################################################################################

def add_particle(i, j):
    grid[i,j] += 1
    check_avalanche_condition(i, j)

###########################################################################################################

def remove_particle(i, j):
    grid[i,j] -= 1
    check_avalanche_condition(i, j)

###########################################################################################################


def hop():
    # this is the ugliest but fastest way
    while True:
        i_start, j_start = np.random.randint(0, width), np.random.randint(0, height)
        # start does not lie in shadow and there actually is sand there
        if not is_in_shadow(i_start, j_start) and grid[i_start, j_start] > 0:
            break

    # particle is in the air, unitl it finds a cell to be stored

    while True:
        # move by wind
        
        this_wind_x = wind_x
        i_end, j_end = (i_start + this_wind_x)%width, (j_start + wind_y)%height

        # particle lands in shadow: 100% is deposited: extension to werner
        if is_in_shadow(i_end, j_end):
            break
        else:
            # propability depents whether new position is sand (>0) or surface (0)
            prob = p_surface if grid[i_end, j_end] == 0 else p_sand
            # set down at new position?
            if np.random.rand() < prob:
                break
    
    # move the particle to the valid position
    remove_particle(i_start, j_start)
    add_particle(i_end, j_end)

###########################################################################################################
        
def is_in_shadow(i, j):
    delta_x = np.arange(1, width)  
    left_indizes = (i - delta_x) % width 
    delta_height = grid[left_indizes, j] - grid[i, j]  
    return np.any(delta_height >= delta_x*shadow_coeff)



###########################################################################################################
###########################################################################################################
###########################################################################################################

def main():
    fig, ax = plt.subplots()
    im = ax.imshow(grid, cmap='YlOrBr', vmin=0, vmax=10)
    plt.colorbar(im)


    def update(frame):
        for _ in range(1000):
            hop()
        im.set_array(grid)
        return [im]

    # Create the animation
    anim = FuncAnimation(fig, update, frames=300, interval=40, blit=True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()