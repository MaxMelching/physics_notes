import numpy as np
import math as m


def gamma(v: float) -> float:
    # return 1./np.sqrt(1. - v**2)
    return 1./m.sqrt(1. - v**2)


def doppler(v: float) -> float:
    # return m.sqrt((1. + v)/(1. - v))
    return (1. + v)*gamma(v)


def add_velocities(v1: float, v2: float) -> float:
    # return v1 + v2
    return (v1 + v2) / (1. + v1 * v2)


def get_half_velocity(v: float) -> float:
    # https://www.wolframalpha.com/input?i=solve+2*v%2F%281%2Bv%5E2%29%3D%3Dx+for+v
    # return (1. + np.sqrt(1 - v**2))/v, (1. - np.sqrt(1 - v**2))/v  # Testing solutions
    # return (1.-np.sqrt(1-v**2))/v
    return (1.-m.sqrt(1-v**2))/v
    # TODO: make distinction with input smaller/greater than zero?


def lorentz(t: float, x: float, v: float) -> tuple[float, float]:
    # return gamma(v)*(t - v*x), gamma(v)*(x - v*t)
    # return gamma(v)*(x - v*t), gamma(v)*(t - v*x)
    return -gamma(v)*(x - v*t), gamma(v)*(t - v*x)  # Minus for intended direction


def draw_clocks(
    v1=0.,
    v2=0.,
    tstart=1.,
    tend=2.,
    t0=1.2,
    t0prime=1.2,
    init_sep=1.,
    referee=True,
) -> str:
    # def obs1(t, x):
    #     return lorentz(t, x, v1)

    # def obs2(t, x):
    #     return lorentz(t, x, v2)

    # if referee:
    #     v_ref = get_half_velocity(add_velocities(v1, v2))

    #     def ref(t, x):
    #         return lorentz(t, x, v_ref)
        
    
    def obs1(t):
        return lorentz(t, init_sep/2., v1)

    def obs2(t):
        return lorentz(t, -init_sep/2., v2)

    if v1 != -v2:
        v_ref = get_half_velocity(add_velocities(v1, v2))
    else:
        v_ref = 0.
    def ref(t):
        return lorentz(t, 0, v_ref)


    def light_l_to_r(t, t0):
        # return (-init_sep/2. + v1 * t0 + t, t0 + t)
        obs_pos = obs1(t0)
        return obs_pos[0] + t, obs_pos[1] + t
        # obs_pos = obs1(t0)
        # travel = lorentz(obs_pos[1] + t, obs_pos[0], v1)
        # return travel[0], travel[1]

    def light_r_to_l(t, t0):
        # return (init_sep/2. + v2 * t0 - t, t0 + t)
        obs_pos = obs2(t0)
        return obs_pos[0] - t, obs_pos[1] + t
        # obs_pos = obs2(t0)
        # travel = lorentz(obs_pos[1] - t, obs_pos[0], v2)
        # return travel[0], travel[1]
    
    # -- No Lorentz transform needed, light trajectory is always the same



    # This here seems to work
    # t0, t0prime = t0prime, t0  # Testing
    tau = (2 + v2 * t0 - v1 * t0) / (1 - v2)  # t0 = t_-
    tauprime = (2 - v1 * t0prime + v2 * t0prime) / (1 + v1)
    # tau, tauprime = tauprime, tau  # Not needed if we switch t0 and t'0
    t_return = (2 - v1 * (t0 + tau) + v2 * (t0 + tau)) / (1 + v1)  # = t_+; Not sure why +v2 * t0prime, but works
    tprime_return = (2 + v2 * (t0prime + tauprime) - v1 * (t0prime + tauprime)) / (1 - v2)  # = t'_+
    # t_return, tprime_return = tprime_return, t_return  # Not needed if we switch t0 and t'0


    # v_rel = add_velocities(v1, -v2)
    v_rel = add_velocities(-v1, v2)
    # v_rel = add_velocities(v1, v2)
    k = doppler(v_rel)

    tau = k*t0prime
    tauprime = k*t0
    # tau = t0 + k*t0
    # tauprime = t0prime + k*t0prime
    
    t_return = k*tauprime
    tprime_return = k*tau
    # t_return = tauprime + k**2*t0
    # tprime_return = tau + k**2*t0prime
    # t_return = tauprime + k*tau
    # tprime_return = tau + k*tauprime

    t = (t0 + t_return) / 2
    tprime = (t0prime + tprime_return) / 2
    # Hmmm, but drawing them on line makes them not lie exactly between t_+, t_-
    # -> just take average of the t_+, t_- coordinates?
    # -> problem: then lines of simultaneity are not parallel anymore...
    
    # t_coords = tuple((np.array(light_l_to_r(0, t0)) + np.array(light_r_to_l(t_return, t0 + tau))) / 2)
    # tprime_coords = tuple((np.array(light_r_to_l(0, t0prime)) + np.array(light_l_to_r(tprime_return, t0prime + tauprime))) / 2)


    t_coords = (light_l_to_r(0, t0)[0] + light_r_to_l(t_return, t0 + tau)[0]) / 2, (light_l_to_r(0, t0)[1] + light_r_to_l(t_return, t0 + tau)[1]) / 2
    tprime_coords = (light_r_to_l(0, t0prime)[0] + light_l_to_r(tprime_return, t0prime + tauprime)[0]) / 2, (light_r_to_l(0, t0prime)[1] + light_l_to_r(tprime_return, t0prime + tauprime)[1]) / 2
    
    

#     return fr'''
# \begin{{tikzpicture}}
#     % Draw lines of observers
#     \draw[->, thick] {obs1(tstart)} -- {obs1(tend)};
#     \draw[->, thick] {obs2(tstart)} -- {obs2(tend)};
#     \draw[->, thick] {ref(tstart)} -- {ref(tend)};

#     % Draw trajectories of light
#     \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(tau, t0)};
#     % \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(0.5*t0*(1. + doppler(v_rel)**2), t0)}; % Light travel time
#     % \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(gamma(v2)*tau, t0)}; % Light travel time
#     % \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r((tau - t0)/gamma(v_rel), t0)}; % Light travel time
#     \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0 + tau)} -- {light_r_to_l(t_return, t0 + tau)};
#     \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tauprime, t0prime)};
#     \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0prime + tauprime)} -- {light_l_to_r(tprime_return, t0prime + tauprime)};

#     % Draw lines of simultaneity
#     \draw[thick, blue] {obs1(t)} -- {obs2(tprime)};
#     %\draw[thick, blue] {t_coords} -- {tprime_coords};
#     \draw[thick, blue] {obs1(t0prime + tauprime)} -- {obs2(t0 + tau)};

#     % Make labels
#     \draw {obs1(tstart)} node[left] {{$\mathcal{{O}}$}};
#     \draw {obs1(t0)} node[left] {{$t_-$}};
#     \draw {obs1(t0 + tau + t_return)} node[left] {{$t_+$}};
#     \draw {obs2(t0 + tau)} node[right] {{$\tau'$}};
#     \draw {obs2(tprime)} node[right] {{$t'$}};
#     %\draw {tprime_coords} node[right] {{$t'$}};
#     \draw {obs2(tstart)} node[right] {{$\mathcal{{O}}'$}};
#     \draw {obs2(t0prime)} node[right] {{$t'_-$}};
#     \draw {obs2(t0prime + tauprime + tprime_return)} node[right] {{$t'_+$}};
#     \draw {obs1(t0prime + tauprime)} node[left] {{$\tau$}};
#     \draw {obs1(t)} node[left] {{$t$}};
#     %\draw {t_coords} node[left] {{$t$}};
#     \draw {ref(tstart)} node[below] {{$\mathcal{{R}}$}};
# \end{{tikzpicture}}
#     '''
    return fr'''
\begin{{tikzpicture}}
    % Draw lines of observers
    \draw[->, thick] {obs1(tstart)} -- {obs1(tend)};
    \draw[->, thick] {obs2(tstart)} -- {obs2(tend)};
    \draw[->, thick] {ref(tstart)} -- {ref(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {obs1(t0)} -- {obs2(t0 + tauprime)};
    \draw[->, thick, black!10!yellow] {obs2(t0 + tauprime)} -- {obs1(t0 + tauprime + t_return)};
    \draw[->, thick, black!10!yellow] {obs2(t0prime)} -- {obs1(t0prime + tau)};
    \draw[->, thick, black!10!yellow] {obs1(t0prime + tau)} -- {obs2(t0prime + tau + tprime_return)};

    % Draw lines of simultaneity
    \draw[thick, blue] {obs1(t + (t0 + tauprime)/2)} -- {obs2(tprime + (t0prime + tau)/2)};
    \draw[thick, blue] {obs1(t0 + tauprime)} -- {obs2(t0prime + tau)};

    % Make labels
    \draw {obs1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {obs1(t0)} node[left] {{$t_-$}};
    \draw {obs1(t0 + tau + t_return)} node[left] {{$t_+$}};
    \draw {obs2(t0 + tau)} node[right] {{$\tau'$}};
    \draw {obs2(tprime + (t0prime + tau)/2)} node[right] {{$t'$}};
    \draw {obs2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {obs2(t0prime)} node[right] {{$t'_-$}};
    \draw {obs2(t0prime + tauprime + tprime_return)} node[right] {{$t'_+$}};
    \draw {obs1(t0prime + tauprime)} node[left] {{$\tau$}};
    \draw {obs1(t + (t0 + tauprime)/2)} node[left] {{$t$}};
    \draw {ref(tstart)} node[below] {{$\mathcal{{R}}$}};
\end{{tikzpicture}}
    '''
#     return fr'''
# \begin{{tikzpicture}}
#     % Draw lines of observers
#     \draw[->, thick] {obs1(tstart)} -- {obs1(tend)};
#     \draw[->, thick] {obs2(tstart)} -- {obs2(tend)};
#     \draw[->, thick] {ref(tstart)} -- {ref(tend)};

#     % Draw trajectories of light
#     \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(tau, t0)};
#     \draw[->, thick, black!10!yellow] {light_r_to_l(0, tau)} -- {light_r_to_l(t_return, tau)};
#     \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tauprime, t0prime)};
#     \draw[->, thick, black!10!yellow] {light_l_to_r(0, tauprime)} -- {light_l_to_r(tprime_return, tauprime)};

#     % Draw lines of simultaneity
#     \draw[thick, blue] {obs1(t)} -- {obs2(tprime)};
#     %\draw[thick, blue] {t_coords} -- {tprime_coords};
#     \draw[thick, blue] {obs1(t0prime + tauprime)} -- {obs2(t0 + tau)};

#     % Make labels
#     \draw {obs1(tstart)} node[left] {{$\mathcal{{O}}$}};
#     \draw {obs1(t0)} node[left] {{$t_-$}};
#     \draw {obs1(t0 + tau + t_return)} node[left] {{$t_+$}};
#     \draw {obs2(t0 + tau)} node[right] {{$\tau'$}};
#     \draw {obs2(tprime)} node[right] {{$t'$}};
#     %\draw {tprime_coords} node[right] {{$t'$}};
#     \draw {obs2(tstart)} node[right] {{$\mathcal{{O}}'$}};
#     \draw {obs2(t0prime)} node[right] {{$t'_-$}};
#     \draw {obs2(t0prime + tauprime + tprime_return)} node[right] {{$t'_+$}};
#     \draw {obs1(t0prime + tauprime)} node[left] {{$\tau$}};
#     \draw {obs1(t)} node[left] {{$t$}};
#     %\draw {t_coords} node[left] {{$t$}};
#     \draw {ref(tstart)} node[below] {{$\mathcal{{R}}$}};
# \end{{tikzpicture}}
#     '''


def to_latex_file(s: str) -> None:
    ...


if __name__ == '__main__':
    print(get_half_velocity(0.5))
    
    print(draw_clocks(
        tend=5,
    ))
    
    print(draw_clocks(
        v1=0.3,
        v2=0.5,
        tend=5,
    ))
    
    print(draw_clocks(
        v1=-0.3,
        v2=0.3,
        tend=5,
    ))

    print(draw_clocks(
        v1=-0.1,
        v2=0.4,
        tstart=-1,
        tend=7,
        # tstart=0,
        # tend=8,
        # init_sep=0,
        # t0=2,
        # t0prime=2,
        # t0=-0.5,
        # t0prime=-0.5,
    ))