#
# Created on:  Tue Oct 08 2024
# By:  Lukas Mettler
# https://lemettler.github.io
#
# Adaptation of the Bishop et al. (2002) Desert Dune Model. 
# Wind direction is normal distributed.
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from tqdm import tqdm
from datetime import datetime
import sys, os, argparse
sys.setrecursionlimit(10_000)
timestamp = datetime.now().strftime("%y-%m-%d-%H:%M:%S")


# initial values
num_frames = 300
steps_per_frame = 500
grid_size = 100
uniform_initialization = True
H = 5
delta_H = 3
delta_avalanche = 3
shadow_angle = 15
l_coeff0, l_coeff1, l_coeff2 = 5, 0.4, 0.002
wind_std = 0 #deg
p_sand = 0.6
p_surface = 0.4
p_shadow = 1 
p_shadow_erosion = 0 
file_path = ''
restore_file = False


# declaration
width, height = None, None
shadow_coeff = None
h_avg = None
grid =  None


def check_avalanche_condition(i, j):
    '''
    Test the neighborhood of (i,j) for the erosion/avalanche condition, that the height difference exeeds a limit delta_avalanche.
    For this implementation we're using a deterministic Moore neigborhood, that also includes the offdiagonal elements.
    The delta_avalanche is approximated to be the same as for the Neumann neighbors and the direction of the max slope is chosen.
    If the avalanche condition is met, the particles are moved.

    - Shadow: particles in shadow can have special erosion probability p_shadow_erosion.
    ***
    - i, j: current position
    '''

    if is_in_shadow(i,j):
        if np.random.rand() < p_shadow_erosion:
            return None
        
    # get the neigboring indizes with periodic boudary condition (modulo)
    l = (i+1)%width
    r = (i-1)%width
    u = (j+1)%height
    d = (j-1)%height

    #bordering_cell_indizes = (l,j), (r,j), (i,u), (i,d) # neumann
    bordering_cell_indizes = (l,j), (r,j), (i,u), (i,d), (l, u), (r, u), (r, d), (l, d) #moore
    
    # get the steepest slope direction
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

def remove_particle(i, j):
    grid[i,j] -= 1
    check_avalanche_condition(i, j)

###########################################################################################################

def calculate_hop_length(i, j):
    '''
    Length a particle gets transported. From Bishop (2002) Eq. 1.
    - Linear term: windward slope and aysmmetric profiles
    - Quadratic term: dune height
    ***
    - i, j: current erosion/last hop position
    ***
    - return: length (float! -> convert to int before applying to grid!)
    '''
    
    h_ref = h_avg - 0.5*np.mean(np.abs(grid - h_avg))
    this_l = l_coeff0  + l_coeff1 * (grid[i, j] - h_ref)
    
    if grid[i, j] >= h_ref:
        this_l += l_coeff2 * (grid[i, j] - h_ref)**2

    return this_l

###########################################################################################################

def pickup_sand():
    '''
    The by far ugliest and at the same time fastest way to do this!
    ***
    - return: position that is not in shadow and contains sand. (int, int)
    '''
    while True:
        i_start, j_start = np.random.randint(0, width), np.random.randint(0, height)
        # start does not lie in shadow and there actually is sand there
        if not is_in_shadow(i_start, j_start) and grid[i_start, j_start] > 0:
            break
    return i_start, j_start

###########################################################################################################

def get_hop_end(i_hop_start, j_hop_start, wind_angle):
    '''
    Position the particle gets transported to, depending on start position and angle.
    ***
    - i_hop_start, j_hop_start: start position of hop. Site of erosion/last hop end.
    - wind_angle: direction the particle gets transported to.
    ***
    - return: end position
    '''
    
    this_hop_length = calculate_hop_length(i_hop_start, j_hop_start)
    i_end = int(i_hop_start + np.cos(wind_angle)*this_hop_length) % width
    j_end = int(j_hop_start + np.sin(wind_angle)*this_hop_length) % height
    return i_end, j_end

###########################################################################################################

def hop():
    '''
    Move one particle until it finds a position to set down and all avalanches are done.
    '''

    i_start, j_start = pickup_sand()
    i_hop_start, j_hop_start = i_start, j_start

    # particle is in the air, unitl it finds a cell to be stored
    remove_particle(i_start, j_start)

    # mean wind direction is theta = 0, along i axis; wind directions is the same for one hop
    wind_angle = np.deg2rad(wind_std)*np.random.randn()

    while True:
        # move by wind
        i_end, j_end = get_hop_end(i_hop_start, j_hop_start, wind_angle)
        i_hop_start, j_hop_start = i_end, j_end

        # propability depends whether new position is sand (>0) or surface (0)
        prob = p_surface if grid[i_end, j_end] == 0 else p_sand
        prob = p_shadow if is_in_shadow(i_end, j_end) else prob
        # set down at new position?
        if np.random.rand() < prob:
            break 
    
    # move the particle to the valid position
    add_particle(i_end, j_end)

###########################################################################################################
        
def is_in_shadow(i, j):
    '''
    Test if a particle is protected from the wind by a particle to the left (mean wind direction).
    IMPORTANT: DIFFERENT WIND DIRECTIONS ARE NOT ACCOUNTED FOR!
    ***
    - i, j: current position
    ***
    - return: wether particle is in shadow or not (bool)
    '''
    delta_x = np.arange(1, width)  
    left_indizes = (i - delta_x) % width 
    delta_height = grid[left_indizes, j] - grid[i, j]  
    return np.any(delta_height >= delta_x*shadow_coeff)

###########################################################################################################

def init_grid():
    '''
    Initialize the grid (global parameter).
    '''
    global width, height, shadow_coeff, h_avg, grid, uniform_initialization, file_path

    width, height = grid_size, grid_size 

    if uniform_initialization: 
        grid = H*np.ones((height, width))  
    else:
        grid = np.random.randint(max(0,H-delta_H), H+delta_H+1, size=(width, height))


    if restore_file:
        file_path = input('File path: ')
        try:
            grid = np.loadtxt(file_path)
            width, height = grid.shape
        except:
            print('Could not load file. Continuing as without...')


    h_avg = np.mean(grid)
    shadow_coeff = np.tan(np.deg2rad(shadow_angle))


###########################################################################################################


def parse_arguments():
    global num_frames, steps_per_frame, grid_size, uniform_initialization, H, delta_H
    global delta_avalanche, shadow_angle, l_coeff0, l_coeff1, l_coeff2, wind_std
    global p_sand, p_surface, p_shadow, p_shadow_erosion, restore_file


    parser = argparse.ArgumentParser(description="Bishop2002 sand dune CA with wind.")
    
    parser.add_argument("--num-frames",                 type=int,       default=num_frames,             help="Number of frames")
    parser.add_argument("--steps-per-frame",            type=int,       default=steps_per_frame,        help="Steps per frame")
    parser.add_argument("--grid-size",                  type=int,       default=grid_size,              help="Grid size")
    parser.add_argument("--uniform-initialization",     type=bool,      default=uniform_initialization, help="Uniform initialization: True / Randomized: False")
    parser.add_argument("--H",                          type=float,     default=H,                      help="Average intialization height")
    parser.add_argument("--delta-H",                    type=float,     default=delta_H,                help="Max Deviation for randomized intialization")
    parser.add_argument("--delta-avalanche",            type=float,     default=delta_avalanche,        help="Height difference that causes an avalanche")
    parser.add_argument("--shadow-angle",               type=float,     default=shadow_angle,           help="Shadow angle")
    parser.add_argument("--l-coeff0",                   type=float,     default=l_coeff0,               help="Hop coeff: const ~ min length")
    parser.add_argument("--l-coeff1",                   type=float,     default=l_coeff1,               help="Hop coeff: linear ~ dune slope")
    parser.add_argument("--l-coeff2",                   type=float,     default=l_coeff2,               help="Hop coeff: quadratic ~ dune height")
    parser.add_argument("--wind-std",                   type=float,     default=wind_std,               help="Wind standard deviation (deg)")
    parser.add_argument("--p-sand",                     type=float,     default=p_sand,                 help="Probability of hopping sand being stored on sand")
    parser.add_argument("--p-surface",                  type=float,     default=p_surface,              help="Probability of hopping sand being stored on stone")
    parser.add_argument("--p-shadow",                   type=float,     default=p_shadow,               help="Probability of hopping sand being stored in shadow (Overrides sand and stone prob.)")
    parser.add_argument("--p-shadow-erosion",           type=float,     default=p_shadow_erosion,       help="Probability of erosion/avalanches in shadow")
    parser.add_argument("--restore-file",               type=bool,      default=restore_file,           help="Restore a previous state")

    args = parser.parse_args()

    # Update global variables with parsed arguments
    num_frames = args.num_frames
    steps_per_frame = args.steps_per_frame
    grid_size = args.grid_size
    uniform_initialization = args.uniform_initialization
    H = args.H
    delta_H = args.delta_H
    delta_avalanche = args.delta_avalanche
    shadow_angle = args.shadow_angle
    l_coeff0, l_coeff1, l_coeff2 = args.l_coeff0, args.l_coeff1, args.l_coeff2
    wind_std = args.wind_std
    p_sand = args.p_sand
    p_surface = args.p_surface
    p_shadow = args.p_shadow
    p_shadow_erosion = args.p_shadow_erosion
    restore_file = args.restore_file



###########################################################################################################
###########################################################################################################    
###########################################################################################################

def main():

    parse_arguments()
    init_grid()

    title = f'{timestamp}' + '\n\n' 
    title += f'Wind: $\sigma_W = {wind_std}°$' + '\n'
    title += f'Hop coefficients: $L_{{0,1,2}}={l_coeff0},{l_coeff1},{l_coeff2}$' + '\n'
    title += f'Surfaces: $P_{{sand}} = {p_sand}$; $P_{{stone}} = {p_surface}$' + '\n'
    title += f'Shadow: $\Theta = {shadow_angle}°$; $P_{{shadow}} = {p_shadow}$; $P_{{erosion}} = {p_shadow_erosion}$' + '\n'
    title += f'Sand: $<h>={h_avg}$; $\Delta h = {delta_avalanche}$' + '\n'
    title += f'Grid: $(w,h) = ({width}, {height})$' + '\n'
    title += f'Anim: frames $= {num_frames}$; steps/frame $= {steps_per_frame}$'

    # plot
    fig, ax = plt.subplots(figsize=(11, 8))
    #im = ax.imshow(grid, cmap='YlOrBr')
    im = ax.imshow(grid, cmap='YlOrBr', vmin=max(0, H - 3*delta_H), vmax=H*2.5)

    #cbar = plt.colorbar(im, location='left', shrink=0.4)
    #cbar.set_label('Sand height')
    plt.xticks([])
    plt.yticks([])
    plt.text(1.01*width, 0.6*height, title, fontsize=9)
    plt.subplots_adjust(left=0.65)
    plt.tight_layout()

    # output directory and file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, 'out')
    data_dir = os.path.join(output_dir, 'data')
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{timestamp}.gif')
    data_file = os.path.join(data_dir, f'{timestamp}.gz')

    # called for each frame
    def update(frame):
        for _ in range(steps_per_frame):
            hop()
        
        im.set_array(grid)
        #im.autoscale()
        #im.set_clim(min(max(0,H-delta_H*2.5), grid.min()), vmax=max(grid.max(), 2.5*H))
        return [im]

    # progress bar
    with tqdm(total=num_frames, desc="Calculating frames", unit="frame") as pbar:
        def update_with_progress(frame):
            result = update(frame) 
            pbar.update(1)
            return result

        anim = FuncAnimation(fig, update_with_progress, frames=num_frames, interval=40, blit=True)
        # save animation
        writer = FFMpegWriter(fps=25, metadata=dict(artist='Me'), bitrate=1800)
        anim.save(output_file, writer=writer)
        print()
        print(f'Written to: {output_file}!')

        np.savetxt(data_file, grid)
        print(f'Written to: {data_file}!')


if __name__ == "__main__":
    main()
