import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys 
sys.setrecursionlimit(10_000)


'''
New
+++++++++++++++++++++++
- Particle hops into shadow: storage with probability p_shadow = 1

+ moore neighborhoods (include offdiagonal corners)
 
- wind speedup like Bishop (2002); very close to Momiji (2000)
'''


# Set up the simulation parameters
H = 5
width, height = 100, 100 
#grid = H*np.ones((height, width))  
grid = np.random.randint(3, 8, size=(width, height))
delta_avalanche = 3


#hop_length =  1
l_coeff0, l_coeff1, l_coeff2 = 5, 0.4, 0.002

p_sand = 0.6
p_surface = 0.4
p_shadow = 1 
p_shadow_erosion = 0 #0.1 #probability to check avalanche condition for particle in shadow

shadow_coeff = np.tan(np.deg2rad(15))
h_avg = np.mean(grid)


###########################################################################################################

def check_avalanche_condition(i, j):

    if is_in_shadow(i,j):
        if np.random.rand() < p_shadow_erosion:
            return None
        
    l = (i+1)%width
    r = (i-1)%width
    u = (j+1)%height
    d = (j-1)%height
    #bordering_cell_indizes = (l,j), (r,j), (i,u), (i,d)
    bordering_cell_indizes = (l,j), (r,j), (i,u), (i,d), (l, u), (r, u), (r, d), (l, d)


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

def calculate_hop_length(i, j):
    # i, j of current hop or erosion site
    # linear: windward slope surface and asymetric profiles
    # nonlinear: dune-height
    # based on Bishop (2002) eq 1

    h_ref = h_avg - 0.5*np.mean(np.abs(grid - h_avg))
    this_l = l_coeff0  + l_coeff1 * (grid[i, j] - h_ref)
    if grid[i, j] >= h_ref:
        this_l += l_coeff2 * (grid[i, j] - h_ref)**2

    return int(this_l)

###########################################################################################################

def hop():
    # this is the ugliest but fastest way
    while True:
        i_start, j_start = np.random.randint(0, width), np.random.randint(0, height)
        # start does not lie in shadow and there actually is sand there
        if not is_in_shadow(i_start, j_start) and grid[i_start, j_start] > 0:
            break

    # particle is in the air, unitl it finds a cell to be stored
    remove_particle(i_start, j_start)

    i_hop_start, j_hop_start = i_start, j_start
    while True:
        # move by wind
        
        this_hop_length = calculate_hop_length(i_hop_start, j_hop_start)

        i_end, j_end = (i_hop_start + this_hop_length) % width, j_hop_start

        # propability depents whether new position is sand (>0) or surface (0)
        prob = p_surface if grid[i_end, j_end] == 0 else p_sand
        prob = p_shadow if is_in_shadow(i_end, j_end) else prob
        # set down at new position?
        if np.random.rand() < prob:
            break 
        else:
            # next hop starts at current cell
            i_hop_start, j_hop_start = i_end, j_end
    
    # move the particle to the valid position
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
