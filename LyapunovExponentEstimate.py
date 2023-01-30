import numpy as np
import matplotlib.pyplot as plt
from scipy import fft, signal, integrate

from mpl_toolkits import mplot3d
from tqdm import tqdm

#contains a variety of methods for analyzing chaotic behavior
#solve for Lyapunov exponent for experimental data
#find dominant periods in data based on FFT
#solve Lorenz and Roessler attractors

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    cos_this_angle = np.clip(np.dot(unit_vector(v1), unit_vector(v2)), -1, 1)
    return np.arccos(cos_this_angle)

def find_periods(x, step_size, plot=False, max_periods=10):
    '''
    Finds period based on FFT periodogram
    args: x = 1D array with data
          step_size = timestep btwn data points
    '''
    
    my_signal = np.copy(x)
    my_signal = np.delete(my_signal, slice(0, int(my_signal.size/3.)))
    my_signal -= np.mean(my_signal)

    f, Pxx = signal.periodogram(my_signal, fs = 1./step_size, window='hanning', scaling='spectrum')

    my_periods = []
    for amp_arg in np.argsort(np.abs(Pxx))[::-1][0:max_periods]:

        period = 1 / f[amp_arg]
        no_similar_period = True

        if np.abs(Pxx[amp_arg]) > 0.3:
            for i in my_periods:
                if np.abs(period-i) < 0.1:
                    no_similar_period = False

            if no_similar_period:
                my_periods.append(period)

    if plot:
        run_duration = my_signal.size*step_size
        # plot signal over 5 of most common period
        fig = plt.figure(figsize = (7, 5))
        plt.plot(np.arange(0, run_duration, step_size), my_signal)
        
        if len(my_periods) != 0 and my_periods[0]*5< run_duration and my_periods[0] != 0:
            plt.xlim(0, 5*my_periods[0])
        else:
            plt.xlim(0, run_duration)

        #plot frequency domain spectrum
        plt.figure(figsize = (7,5))
        plt.plot(f, Pxx)
        plt.xlim(0, 1)
        plt.xlabel('Frequency (cycles/t.u.)')
        plt.ylabel('Spectrum Amplitude')
    
    return my_periods

def lyapunov_solve(x, step_sz, period, evolv_t_orbits=1, scalmx=0.5, scalmn=0.05, 
                   anglmx=0.2, print_debug=0, plot_replace=0, plot_convergence=0):
    '''
    Estimate the largest Lyapunov exponent of input experimental data of dimension
    at least 3. If the data settles on a stable attracting point, break calc and
    return -100. If the data is periodic for a long time, break calc and return 0.

    Args: x = the experimental data / trajectory, dim >=3
          step_sz = time btwn datapts in trajectory
          period = period or indicative timescale of data
          evolv_t_orbits = # of orbits between vector replacements
          scalmx = max length of vector replacement
          scalmn = min length of vector replacement
          anglmx = max angle of vector replacement
          print_debug = print vals for debugging
          plot_replace = plot the vector replacement process for data dim=3
          plot_convergence = plot convergence of lambda_1 with time for data dim=3
    '''

    #system, calc. params
    npt = x.shape[0] # num samples
    orbit = int(period/step_sz) # samples within an orbit
    fid_pt = x[0] # I.C. of orbit
    evolv_t_steps = int(evolv_t_orbits*orbit) #steps btwn replacements

    if print_debug:
       print("npt ", npt, " orbit ", orbit, " evolv_t_steps ", evolv_t_steps)

    #first, find the closest pt.
    closest_pt_ind = orbit
    closest_pt_dist = np.linalg.norm(fid_pt - x[orbit])

    for i in range(orbit+1, npt-evolv_t_steps-1):
        this_dist = np.linalg.norm(fid_pt - x[i])
        if scalmn < this_dist < closest_pt_dist:
            closest_pt_dist = this_dist
            closest_pt_ind = i

    if closest_pt_dist < scalmn:
        return -100 #the datapts are too close, stable attracting point

    #calculation, plotting params
    closest_pt = x[orbit]
    lambda_1 = 0
    i_max = npt-evolv_t_steps
    no_div_steps = 0
    if plot_replace:
        fid_pt_plt = []
        old_l_plt = []
        new_l_plt = []
    if plot_convergence:
        lambda_1_evolv = []
        t = []

    #now, travel along the fiducial trajectory, measuring divergence and
    #performing vector replacement every evolv_t_steps
    for i in tqdm(range(evolv_t_steps, i_max, evolv_t_steps)):

        if print_debug:
            print("this ind ", i, " closest pt ind ", closest_pt_ind)

        # calculate current l
        l = np.linalg.norm(closest_pt - fid_pt) #L(t_{k-1})
        
        # increment the trajectories by evolv_t_steps
        # i < i_max = npt-evolv_t_steps to avoid inde
        fid_pt = x[i]
        closest_pt_ind += evolv_t_steps
        closest_pt = x[closest_pt_ind]

        # update lyapunov calc sum
        l_prime = closest_pt - fid_pt
        l_prime_dist = np.linalg.norm(closest_pt - fid_pt) #L'(t_k)
        if l == 0:
          return -100
        else:
          lambda_1 += np.log2(l_prime_dist/l)
        this_t = i*step_sz
        this_lambda_1 = lambda_1/(this_t)

        # check for consecutive lack of divergence, break early if so
        if l_prime_dist == l:
            no_div_steps += 1
        else:
            no_div_steps = 0
        if no_div_steps == 10:
            return this_lambda_1

        # plot and print debug code
        if print_debug:
            print("i ", i)
            print("l ", l, " l' ", l_prime)
            print("lambda1 ", lambda_1)
        if plot_convergence:
            lambda_1_evolv.append(this_lambda_1)
            t.append(this_t)
        if plot_replace:
            l_prime_pt = closest_pt

        #look for replacement pt
        scalmx_mult = 1.0
        closest_l_angle = anglmx
        search = 1

        #loop thru pts on trajectory that are at least an orbit away from
        #fiducial point.
        #if we don't find a neighbor pt, loosen our requirements
        while (search):
            # loop thru pts
            for j in range(npt-evolv_t_steps-1):
                # don't take points within same orbit of fid pt
                if (i-orbit) < j < (i+orbit):
                    break
                
                this_dist = np.linalg.norm(fid_pt - x[j])

                # choose pt within distance range, and having closest angular dist.
                if scalmn < this_dist < scalmx_mult*scalmx and this_dist < l_prime_dist:
                    this_angle = angle_between(l_prime, x[j]-fid_pt)

                    #found a pt, update our closest_pt values, turn off search
                    if this_angle < closest_l_angle: 
                        closest_l_angle = this_angle
                        closest_pt = x[j]
                        closest_pt_ind = j
                        search = 0

            #loosen requirenents
            if search:
                if scalmx_mult < 5:
                    scalmx_mult +=1
                elif closest_l_angle == anglmx:
                    closest_l_angle = anglmx*2
                else:
                    break #fail to find a closer point, search=1
        
        if plot_replace:
            fid_pt_plt.append(fid_pt);
            old_l_plt.append(l_prime_pt)
            new_l_plt.append(closest_pt)

        # if we didn't find a new point, break if the neighbor is at the end
        # of the trajectory. otherwise, leave as is
        if search:
            if closest_pt_ind + evolv_t_steps >= npt:
                print("no close neighbors, broke calc at", i, "/", npt, " timesteps")
                return lambda_1/(i*step_sz)

    if plot_replace:
        fig = plt.figure(figsize=(20, 15))
        ax = plt.axes(projection='3d')
        ax.set_title('Calculation Replacement Process')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        for i in range(len(fid_pt_plt)):
            this_fid_pt = fid_pt_plt[i]
            this_old_l_pt = old_l_plt[i]
            this_new_l_pt = new_l_plt[i]

            old_l_x = np.asarray([this_fid_pt[0], this_old_l_pt[0]])
            old_l_y = np.asarray([this_fid_pt[1], this_old_l_pt[1]])
            old_l_z = np.asarray([this_fid_pt[2], this_old_l_pt[2]])

            new_l_x = np.asarray([this_fid_pt[0], this_new_l_pt[0]])
            new_l_y = np.asarray([this_fid_pt[1], this_new_l_pt[1]])
            new_l_z = np.asarray([this_fid_pt[2], this_new_l_pt[2]])

            ax.plot(old_l_x, old_l_y, old_l_z, color='royalblue')
            ax.plot(new_l_x, new_l_y, new_l_z, color='coral')

        fid_pt_plt = np.asarray(fid_pt_plt)
        ax.scatter(fid_pt_plt[:,0], fid_pt_plt[:,1], fid_pt_plt[:,2], color='red', marker='o')
        plt.show()

    if plot_convergence:
        fig = plt.figure(figsize=(12,8))
        plt.title('Convergence of Lyapunov Exponent Over Time')
        plt.xlabel('Time')
        plt.ylabel('Lyapunov Exponent')
        plt.plot(t, lambda_1_evolv, color='black')

    return lambda_1/(npt*step_sz)


def lyapunov_solve_unknown(x,step_sz, default_pd=10):
    '''
    Estimates the Lyapunov exponent of experimental
    data. Supplies the x_max and period args to the
    lyapunov_solve function
    args: x = experimental data, dim >=3
          step_sz = timestep btwn datapoints
          default_pd = period if FFT can't find a period
    '''

    period_estimate = find_periods(x[:,0], step_sz, max_periods=5)
    if len(period_estimate) > 0:
        poss_periods = [i for i,v in enumerate(period_estimate) if v>0.5]
        if len(poss_periods) > 0:
            period_estimate = period_estimate[poss_periods[0]]
        else:
            period_estimate = default_pd
    else:
        period_estimate= default_pd

    x_max = 0.1*(np.ptp(x,axis=0)[0])

    return lyapunov_solve(x, step_sz, period_estimate, scalmx=x_max)

# =========================================================================
# example attractors
# =========================================================================

a = 0.15
b = 0.20
c = 10.0

def roessler_eqn(x_3, t0):
  x, y, z = x_3
  return [-y-z, x+a*y, b+z*(x-c)]

def solve_roessler(max_time = 500, step_size = 0.01):

  #choose N random starting points, uniformly distributed [-15, 15)
  np.random.seed(1)
  x0 = -15 + 30*np.random.random(3) 

  #define timesteps [0, mt). 250 timesteps per unit time
  t = np.linspace(0, max_time, int(max_time/step_size))

  #integrate ODE for each random starting point
  x_t = np.asarray(integrate.odeint(roessler_eqn, x0, t))

  return t, x_t

# key params
sigma = 6.5
beta = 8.0/3
rho = 30.0

def lorenz_eqn(x_3, t0):
  x, y, z = x_3
  return [sigma*(y-x), x*(rho-z) - y, x*y - beta*z]

def solve_lorenz(max_time = 500, step_size = 0.01):

  #choose N random starting points, uniformly distributed [-15, 15)
  np.random.seed(1)
  x0 = -15 + 30*np.random.random(3) 

  #define timesteps [0, mt). 250 timesteps per unit time
  t = np.linspace(0, max_time, int(max_time/step_size))

  #integrate ODE for each random starting point
  x_t = np.asarray(integrate.odeint(lorenz_eqn, x0, t))

  return t, x_t
