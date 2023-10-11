import numpy as np



def resting_clocks_v1(tstart, tend, t1):
    """
    Draw clock diagram for two observers at rest.
    """
    def line1(t):
        return (-1, t)
    
    def line2(t):
        return (1, t)
    
    light_travel_time = 2
    t_intersect = t1 + light_travel_time
    tprime_intersect = t1 + light_travel_time
    t_return = t1 + 2 * light_travel_time

    
    return fr'''
\subfloat[Time $t$ measured by observer $\mathcal{{O}}$]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {line1(t1)} -- {line2(t_intersect)};
    \draw[->, thick, black!10!yellow] {line2(t_intersect)} -- {line1(t_return)};

    % Make labels
    \draw (-1, {tstart}) node[left] {{$\mathcal{{O}}$}};
    \draw (0, {(t1 + t_intersect) / 2}) node[below right] {{$t_-$}};
    \draw (0, {(t_intersect + t_return) / 2}) node[above right] {{$t_+$}};
    \draw (1, {t_intersect}) node[right] {{$t$}};
    \draw (1, {tstart}) node[right] {{$\mathcal{{O}}'$}};
\end{{tikzpicture}}
}}\hspace*{{0.1\textwidth}}
%
\subfloat[Time $t'$ measured by observer $\mathcal{{O}}'$]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {line2(t1)} -- {line1(t_intersect)};
    \draw[->, thick, black!10!yellow] {line1(t_intersect)} -- {line2(t_return)};

    % Make labels
    \draw (-1, {tstart}) node[right] {{$\mathcal{{O}}$}};
    \draw (1, {tstart}) node[right] {{$\mathcal{{O}}'$}};
    \draw (0, {(t1 + t_intersect) / 2}) node[below left] {{$t'_-$}};
    \draw (0, {(tprime_intersect + t_return) / 2}) node[above left] {{$t'_+$}};
    \draw (-1, {tprime_intersect}) node[left] {{$t'$}};
\end{{tikzpicture}}
}}\hspace*{{0.1\textwidth}}
%
\subfloat[Comparison of both observers]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {line1(t1)} -- {line2(t_intersect)};
    \draw[->, thick, black!10!yellow] {line2(t_intersect)} -- {line1(t_return)};
    \draw[->, thick, black!10!yellow] {line2(t1)} -- {line1(tprime_intersect)};
    \draw[->, thick, black!10!yellow] {line1(tprime_intersect)} -- {line2(t_return)};

    % Make labels
    \draw (-1, {tstart}) node[left] {{$\mathcal{{O}}$}};
    \draw (1, {t_intersect}) node[right] {{$t$}};
    \draw (1, {tstart}) node[right] {{$\mathcal{{O}}'$}};
    \draw (-1, {tprime_intersect}) node[left] {{$t'$}};
\end{{tikzpicture}}
}}
    '''


def resting_clocks(tstart, tend, t0, x_distance=1.0, x_shift=0.0):
    """
    Draw clock diagram for two observers at rest.

    tstart: starting time for observers
    tend: ending time for observers
    t0: time where light first ray is sent out from both observers
    x_distance: distance of each observer from origin
    x_shift: shift for whole coordinates into positive x-direction
    """
    def line1(t):  # Observer 1
        return (-x_distance + x_shift, t)
    
    def line2(t):  # Observer 2
        return (x_distance + x_shift, t)
    
    def line3(t):  # Referee
        return (0 + x_shift, t)
    
    light_travel_time = 2 * x_distance
    t_intersect = t0 + light_travel_time
    tprime_intersect = t0 + light_travel_time
    t_return = t0 + 2 * light_travel_time

    
    return fr'''
\subfloat[Time $t'$ measured by observer $\mathcal{{O}}$]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {line1(t0)} -- {line2(t_intersect)};
    \draw[->, thick, black!10!yellow] {line2(t_intersect)} -- {line1(t_return)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line1(t0)} node[left] {{$t_-$}};
    \draw {line1(t_return)} node[left] {{$t_+$}};
    \draw {line2(t_intersect)} node[right] {{$t'$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
\end{{tikzpicture}}
}}\hspace*{{0.1\textwidth}}
%
\subfloat[Time $t$ measured by observer $\mathcal{{O}}'$]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {line2(t0)} -- {line1(t_intersect)};
    \draw[->, thick, black!10!yellow] {line1(t_intersect)} -- {line2(t_return)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {line2(t0)} node[right] {{$t'_-$}};
    \draw {line2(t_return)} node[right] {{$t'_+$}};
    \draw {line1(tprime_intersect)} node[left] {{$t$}};
\end{{tikzpicture}}
}}\hspace*{{0.1\textwidth}}
%
\subfloat[Comparison of both observers]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};
    \draw[->, thick] {line3(tstart)} -- {line3(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {line1(t0)} -- {line2(t_intersect)};
    \draw[->, thick, black!10!yellow] {line2(t_intersect)} -- {line1(t_return)};
    \draw[->, thick, black!10!yellow] {line2(t0)} -- {line1(tprime_intersect)};
    \draw[->, thick, black!10!yellow] {line1(tprime_intersect)} -- {line2(t_return)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line2(t_intersect)} node[right] {{$t'$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {line1(tprime_intersect)} node[left] {{$t$}};
    \draw {line3(tstart)} node[below] {{$\mathcal{{R}}$}};
\end{{tikzpicture}}
}}
    '''



def resting_clocks_single_fig(tstart, tend, t0, x_distance=1.0, x_shift=0.0):
    """
    Draw clock diagram for two observers at rest.

    tstart: starting time for observers
    tend: ending time for observers
    t0: time where light first ray is sent out from both observers
    x_distance: distance of each observer from origin
    x_shift: shift for whole coordinates into positive x-direction
    """
    def line1(t):  # Observer 1
        return (-x_distance + x_shift, t)
    
    def line2(t):  # Observer 2
        return (x_distance + x_shift, t)
    
    def line3(t):  # Referee
        return (0 + x_shift, t)
    
    light_travel_time = 2 * x_distance
    t_intersect = t0 + light_travel_time
    tprime_intersect = t0 + light_travel_time
    t_return = t0 + 2 * light_travel_time

    
    return fr'''
\subfloat[Comparison of both observers]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};
    \draw[->, thick] {line3(tstart)} -- {line3(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {line1(t0)} -- {line2(t_intersect)};
    \draw[->, thick, black!10!yellow] {line2(t_intersect)} -- {line1(t_return)};
    \draw[->, thick, black!10!yellow] {line2(t0)} -- {line1(tprime_intersect)};
    \draw[->, thick, black!10!yellow] {line1(tprime_intersect)} -- {line2(t_return)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line1(t0)} node[left] {{$t_-$}};
    \draw {line1(t_return)} node[left] {{$t_+$}};
    \draw {line2(t_intersect)} node[right] {{$t'$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {line2(t0)} node[right] {{$t'_-$}};
    \draw {line2(t_return)} node[right] {{$t'_+$}};
    \draw {line1(tprime_intersect)} node[left] {{$t$}};
    \draw {line3(tstart)} node[below] {{$\mathcal{{R}}$}};
\end{{tikzpicture}}
}}
    '''



def moving_clocks_resting_wrt_each_other(v, tstart, tend, t0, t0prime=None, x_distance=1):
    """
    Draw clock diagram for two moving observers which are at rest with
    respect to each other.

    v: velocity that both observers move with
    tstart: starting time for observers
    tend: ending time for observers
    t0: time where first light ray is sent out from first observer.
    t0prime: time where first light ray is sent out from second observer.
    If it is None, we choose it such that clocks are synchronized.
    x_distance: distance of each observer from origin
    """
    x_distance = np.abs(x_distance)  # Ensure it is positive


    def line1(t):  # Observer 1
        return (-x_distance + v * t, t)
    
    def line2(t):  # Observer 2
        return (x_distance + v * t, t)
    
    def line3(t):  # Referee
        return (v * t, t)

    def light_l_to_r(t, t0):
        return (-x_distance + v * t0 + t, t0 + t)

    def light_r_to_l(t, t0):
        return (x_distance + v * t0 - t, t0 + t)
    

    if t0prime is None:
        # Set value such that clocks are synchronizeds
        t_ref = x_distance / (1 - v)  # Time where light from left observer reaches referee
        tprime_ref = x_distance / (1 + v)  # Time where light from right observer reaches referee
        # t1prime = (2 + v * t1 - t_ref) / v - t_ref
        # t1prime = t1 - t_ref + tprime_ref
        t0prime = t0 + t_ref - tprime_ref  # No idea why this works, found it by testing -> ah, t0 + tref = t0prime + t0ref
    

    # t = (2 - v * t1) / (1 - v)  # t1 = t_-
    t = 2 * x_distance / (1 - v)  # t1 = t_-
    # t_return = (2 + v * (t1 + t)) / (1 + v)  # = t_+
    t_return = 2 * x_distance / (1 + v)  # = t_+
    # tprime = (2 + v * t1) / (1 + v)
    tprime = 2 * x_distance / (1 + v)
    # tprime_return = (2 - v * t1) / (1 - v)
    tprime_return = 2 * x_distance / (1 - v)

    
    return fr'''
\subfloat[Time $t'$ measured by observer $\mathcal{{O}}$]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(t, t0)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0 + t)} -- {light_r_to_l(t_return, t0 + t)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line1(t0)} node[left] {{$t_-$}};
    \draw {line1(t0 + t + t_return)} node[left] {{$t_+$}};
    \draw {line2(t0 + t)} node[right] {{$t'$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
\end{{tikzpicture}}
}}\hspace*{{0.025\textwidth}}
%
\subfloat[Time $t$ measured by observer $\mathcal{{O}}'$]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tprime, t0prime)};
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0prime + tprime)} -- {light_l_to_r(tprime_return, t0prime + tprime)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {line2(t0prime)} node[right] {{$t'_-$}};
    \draw {line2(t0prime + tprime + tprime_return)} node[right] {{$t'_+$}};
    \draw {line1(t0prime + tprime)} node[left] {{$t$}};
\end{{tikzpicture}}
}}\hspace*{{0.025\textwidth}}
%
\subfloat[Comparison of both observers]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};
    \draw[->, thick] {line3(tstart)} -- {line3(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(t, t0)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0 + t)} -- {light_r_to_l(t_return, t0 + t)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tprime, t0prime)};
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0prime + tprime)} -- {light_l_to_r(tprime_return, t0prime + tprime)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line2(t0 + t)} node[right] {{$t'$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {line1(t0prime + tprime)} node[left] {{$t$}};
    \draw {line3(tstart)} node[below] {{$\mathcal{{R}}$}};
\end{{tikzpicture}}
}}
    '''



def moving_clocks_resting_wrt_each_other_single_fig(v, tstart, tend, t0, t0prime=None, x_distance=1):
    """
    Draw clock diagram for two moving observers which are at rest with
    respect to each other.

    v: velocity that both observers move with
    tstart: starting time for observers
    tend: ending time for observers
    t0: time where first light ray is sent out from first observer.
    t0prime: time where first light ray is sent out from second observer.
    If it is None, we choose it such that clocks are synchronized.
    x_distance: distance of each observer from origin
    """
    x_distance = np.abs(x_distance)  # Ensure it is positive


    def line1(t):  # Observer 1
        return (-x_distance + v * t, t)
    
    def line2(t):  # Observer 2
        return (x_distance + v * t, t)
    
    def line3(t):  # Referee
        return (v * t, t)

    def light_l_to_r(t, t0):
        return (-x_distance + v * t0 + t, t0 + t)

    def light_r_to_l(t, t0):
        return (x_distance + v * t0 - t, t0 + t)
    

    if t0prime is None:
        # Set value such that clocks are synchronizeds
        t_ref = x_distance / (1 - v)  # Time where light from left observer reaches referee
        tprime_ref = x_distance / (1 + v)  # Time where light from right observer reaches referee
        # t1prime = (2 + v * t1 - t_ref) / v - t_ref
        # t1prime = t1 - t_ref + tprime_ref
        t0prime = t0 + t_ref - tprime_ref  # No idea why this works, found it by testing -> ah, t0 + tref = t0prime + t0ref
    

    # t = (2 - v * t1) / (1 - v)  # t1 = t_-
    t = 2 * x_distance / (1 - v)  # t1 = t_-
    # t_return = (2 + v * (t1 + t)) / (1 + v)  # = t_+
    t_return = 2 * x_distance / (1 + v)  # = t_+
    # tprime = (2 + v * t1) / (1 + v)
    tprime = 2 * x_distance / (1 + v)
    # tprime_return = (2 - v * t1) / (1 - v)
    tprime_return = 2 * x_distance / (1 - v)

    
    return fr'''
\subfloat[Comparison of both observers]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};
    \draw[->, thick] {line3(tstart)} -- {line3(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(t, t0)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0 + t)} -- {light_r_to_l(t_return, t0 + t)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tprime, t0prime)};
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0prime + tprime)} -- {light_l_to_r(tprime_return, t0prime + tprime)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line1(t0)} node[left] {{$t_-$}};
    \draw {line1(t0 + t + t_return)} node[left] {{$t_+$}};
    \draw {line2(t0 + t)} node[right] {{$t'$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {line2(t0prime)} node[right] {{$t'_-$}};
    \draw {line2(t0prime + tprime + tprime_return)} node[right] {{$t'_+$}};
    \draw {line1(t0prime + tprime)} node[left] {{$t$}};
    \draw {line3(tstart)} node[below] {{$\mathcal{{R}}$}};
\end{{tikzpicture}}
}}
    '''



def resting_clocks_fitted(tstart, tend, t0, v, x_distance=1.0, x_shift=0.0):
    """
    Draw clock diagram for two observers at rest.
    Now with unequal distances and positions

    tstart: starting time for observers
    tend: ending time for observers
    t0: time where light first ray is sent out from both observers
    x_distance: distance of each observer from origin
    x_shift: shift for whole coordinates into positive x-direction
    """
    # Compute times where they shall intersect -> not even needed... We know how to compute this more easily from previous function
    # t_int1, t_int2, t_int3 = (tend - tstart - t0) / 2, (tend - tstart) / 2, (tend - tstart + t0) / 2
    # # Compute shifts. Done by substituting fixed times + locations into Lorentz transform.
    # # Here we do inverse transform, thus we use +v instead of -v
    # gamma = 1 / np.sqrt(1 - v**2)
    # x_shift1, x_shift2, x_shift3 = (x_shift - x_distance + v * t_int1) * gamma, (x_shift + v * t_int2) * gamma, (x_shift + x_distance + v * t_int3) * gamma
    # # Note: here 2 and 3 are switched compared to line names...

    t = 2 * x_distance / (1 - v)
    tprime = 2 * x_distance / (1 + v)
    # Set value such that clocks are synchronizeds
    t_ref = x_distance / (1 - v)  # Time where light from left observer reaches referee
    tprime_ref = x_distance / (1 + v)  # Time where light from right observer reaches referee
    t0prime = t0 + t_ref - tprime_ref  # No idea why this works, found it by testing -> ah, t0 + tref = t0prime + t0ref

    x_shift1, x_shift2, x_shift3 = -x_distance + v * (t0prime + tprime), v * (tend - tstart) / 2, x_distance + v * (t0 + t)

    def line1(t):  # Observer 1
        return (x_shift1, t)
    
    def line2(t):  # Observer 2
        return (x_shift3, t)
    
    def line3(t):  # Referee
        return (x_shift2, t)
    
    light_travel_time = x_shift3 - x_shift1
    t_intersect = t0 + light_travel_time
    tprime_intersect = t0 + light_travel_time
    t_return = t0 + 2 * light_travel_time

    
    return fr'''
\subfloat[Comparison of both observers]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};
    \draw[->, thick] {line3(tstart)} -- {line3(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {line1(t0)} -- {line2(t_intersect)};
    \draw[->, thick, black!10!yellow] {line2(t_intersect)} -- {line1(t_return)};
    \draw[->, thick, black!10!yellow] {line2(t0)} -- {line1(tprime_intersect)};
    \draw[->, thick, black!10!yellow] {line1(tprime_intersect)} -- {line2(t_return)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line2(t_intersect)} node[right] {{$t'$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {line1(tprime_intersect)} node[left] {{$t$}};
    \draw {line3(tstart)} node[below] {{$\mathcal{{R}}$}};
\end{{tikzpicture}}
}}
    '''



def moving_clocks(v1, v2, tstart, tend, t0, t0prime=None):
    """
    Draw clock diagram for two moving observers which are at rest with
    respect to each other.

    If t1prime is None, we choose it such that clocks are synchronized.
    """
    v_avg = (v1 + v2) / 2
    # v_rel = np.abs(v1 - v2)

    def line1(t):  # Observer 1
        return (-1 + v1 * t, t)
        # return tuple([-1, 0] + [v1, 1] / np.sqrt(v1**2 + 1) * t)
    
    def line2(t):  # Observer 2
        return (1 + v2 * t, t)
        # return tuple([1, 0] + [v2, 1] / np.sqrt(v2**2 + 1) * t)
    
    def line3(t):  # Referee
        # if v_avg != 0:
        #     return (v1 * v2 / v_avg * t, t)
        #     # return (t, t * v_avg / v1 / v2)
        # else:
        #     return (t, 0)

        # With formula from Mathematica, input was:
        # Solve[Norm[{v3*t, t}/Norm[{v3*t, t}] - {v2*t, t}/Norm[{v2*t, t}]] == 
        # Norm[{v3*t, t}/Norm[{v3*t, t}] - {v1*t, t}/Norm[{v1*t, t}]], v3, Reals, 
        # Assumptions ->  v1 \[Element] Reals && v2 \[Element] Reals && v3 \[Element] Reals && t \[Element] Reals]
        # if v_avg < 0:
        #     v3 = (-1 + v1 * v2) / (v1 + v2) - np.sqrt(1 + v1**2 + v2**2 + v1**2 * v2**2) / np.abs(v1 + v2)
        #     # print(v3)
        #     return (v3 * t, t)
        # elif v_avg > 0:
        #     v3 = (-1 + v1 * v2) / (v1 + v2) + np.sqrt(1 + v1**2 + v2**2 + v1**2 * v2**2) / np.abs(v1 + v2)
        #     # print(v3)
        #     return (v3 * t, t)
        # else:
        #     return (0, t)
        
        return (v_avg * t, t)

    def light_l_to_r(t, t0):
        return (-1 + v1 * t0 + t, t0 + t)

    def light_r_to_l(t, t0):
        return (1 + v2 * t0 - t, t0 + t)
    

    # v1 /= np.sqrt(v1**2 + 1)
    # v2 /= np.sqrt(v2**2 + 1)
    

    if t0prime is None:
        # Set value such that clocks are synchronizeds
        t_ref = (1 + v_avg * t0 - v1 * t0) / (1 - v_avg)  # Time where light from left observer reaches referee
        # tprime_ref = (1 + v2 * t0 + v1 * t0) / (1 + v2)  # Time where light from right observer reaches referee
        # tprime_ref = 1 / (1 + v2)  # Time where light from right observer reaches referee
        # t0prime = t0 + t_ref - tprime_ref
        t0prime = (t0 + t_ref - 1 / (1 + v_avg)) / (1 + (v2 - v_avg) / (1 + v_avg))
    

    # t = (2 + v2 * t0 - v1 * t0) / (1 - v2)  # t1 = t_-
    # # t = 2 / (1 - v2)  # t1 = t_-
    # # t_return = (2 + v1 * (t0 + t) - v2 * (t0 + t)) / (1 + v2)  # = t_+
    # t_return = (2 + v1 * (t0 + t) + v2 * (t0 + t)) / (1 + v2)  # = t_+; Not sure why +v2 * t0prime, but works
    # # t_return = 2 / (1 + v2)  # = t_+
    # tprime = (2 + v1 * t0prime + v2 * t0prime) / (1 + v2)  # Not sure why +v2 * t0prime, but works
    # # tprime = 2 / (1 + v2)
    # tprime_return = (2 + v2 * (t0prime + tprime) - v1 * (t0prime + tprime)) / (1 - v2)
    # # tprime_return = 2 / (1 - v2)


    # t = (2 + v2 * t0 - v1 * t0) / (1 - v2)  # t1 = t_-
    # tprime = (2 + v1 * t0 + v2 * t0prime) / (1 + v2)  # Not sure why +v2 * t0prime, but works
    # t_return = (2 + v1 * (t0 + t) + v2 * (t0prime + tprime)) / (1 + v2)  # = t_+; Not sure why +v2 * t0prime, but works
    # tprime_return = (2 + v2 * (t0prime + tprime) - v1 * (t0 + t)) / (1 - v2)


    # t = (2 + v_rel * t0) / (1 - v2)  # t1 = t_-
    # # t = 2 / (1 - v2)  # t1 = t_-
    # # t_return = (2 + v1 * (t0 + t) - v2 * (t0 + t)) / (1 + v2)  # = t_+
    # t_return = (2 + 2 * v_avg * (t0 + t)) / (1 + v2)  # = t_+; Not sure why +v2 * t0prime, but works
    # # t_return = 2 / (1 + v2)  # = t_+
    # tprime = (2 + 2 * v_avg * t0prime) / (1 + v2)  # Not sure why +v2 * t0prime, but works
    # # tprime = 2 / (1 + v2)
    # tprime_return = (2 + v_rel * (t0prime + tprime)) / (1 - v2)


    # This here seems to work
    # t0, t0prime = t0prime, t0  # Testing
    tau = (2 + v2 * t0 - v1 * t0) / (1 - v2)  # t0 = t_-
    tauprime = (2 - v1 * t0prime + v2 * t0prime) / (1 + v1)
    # tau, tauprime = tauprime, tau  # Not needed if we switch t0 and t'0
    t_return = (2 - v1 * (t0 + tau) + v2 * (t0 + tau)) / (1 + v1)  # = t_+; Not sure why +v2 * t0prime, but works
    tprime_return = (2 + v2 * (t0prime + tauprime) - v1 * (t0prime + tauprime)) / (1 - v2)  # = t'_+
    # t_return, tprime_return = tprime_return, t_return  # Not needed if we switch t0 and t'0


    t = (t0 + t_return) / 2
    tprime = (t0prime + tprime_return) / 2
    # Hmmm, but drawing them on line makes them not lie exactly between t_+, t_-
    # -> just take average of the t_+, t_- coordinates?
    # -> problem: then lines of simultaneity are not parallel anymore...
    t_coords = tuple((np.array(light_l_to_r(0, t0)) + np.array(light_r_to_l(t_return, t0 + tau))) / 2)
    tprime_coords = tuple((np.array(light_r_to_l(0, t0prime)) + np.array(light_l_to_r(tprime_return, t0prime + tauprime))) / 2)



    
    return fr'''
\subfloat[Time $\tau'$ measured by observer $\mathcal{{O}}$]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(tau, t0)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0 + tau)} -- {light_r_to_l(t_return, t0 + tau)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line1(t0)} node[left] {{$t_-$}};
    \draw {line1(t0 + tau + t_return)} node[left] {{$t_+$}};
    \draw {line2(t0 + tau)} node[right] {{$\tau'$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
\end{{tikzpicture}}
}}\hspace*{{0.025\textwidth}}
%
\subfloat[Time $\tau$ measured by observer $\mathcal{{O}}'$]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tauprime, t0prime)};
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0prime + tauprime)} -- {light_l_to_r(tprime_return, t0prime + tauprime)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {line2(t0prime)} node[right] {{$t'_-$}};
    \draw {line2(t0prime + tauprime + tprime_return)} node[right] {{$t'_+$}};
    \draw {line1(t0prime + tauprime)} node[left] {{$\tau$}};
\end{{tikzpicture}}
}}\hspace*{{0.025\textwidth}}
%
\subfloat[Comparison of both observers]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};
    \draw[->, thick] {line3(tstart)} -- {line3(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(tau, t0)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0 + tau)} -- {light_r_to_l(t_return, t0 + tau)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tauprime, t0prime)};
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0prime + tauprime)} -- {light_l_to_r(tprime_return, t0prime + tauprime)};

    % Draw lines of simultaneity
    \draw[thick, blue] {line1(t)} -- {line2(tprime)};
    %\draw[thick, blue] {t_coords} -- {tprime_coords};
    \draw[thick, blue] {line1(t0prime + tauprime)} -- {line2(t0 + tau)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line2(t0 + tau)} node[right] {{$\tau'$}};
    \draw {line2(tprime)} node[right] {{$t'$}};
    %\draw {tprime_coords} node[right] {{$t'$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {line1(t0prime + tauprime)} node[left] {{$\tau$}};
    \draw {line1(t)} node[left] {{$t$}};
    %\draw {t_coords} node[left] {{$t$}};
    \draw {line3(tstart)} node[below] {{$\mathcal{{R}}$}};
\end{{tikzpicture}}
}}
    '''


    
    # %\draw[->, thick, black!10!yellow] {line1(t0)} -- {line2(tau)};
    # %\draw[->, thick, black!10!yellow] {line2(tau)} -- {line1(t_return)};
    # %\draw[->, thick, black!10!yellow] {line2(t0prime)} -- {line1(tauprime)};
    # %\draw[->, thick, black!10!yellow] {line1(tauprime)} -- {line2(tprime_return)};

    # %\draw[->, thick, black!10!yellow] {line1(t0)} -- {line2(t0prime + tauprime)};
    # %\draw[->, thick, black!10!yellow] {line2(t0prime + tauprime)} -- {line1(t0 + tau + t_return)};
    # %\draw[->, thick, black!10!yellow] {line2(t0prime)} -- {line1(t0 + tau)};
    # %\draw[->, thick, black!10!yellow] {line1(t0 + tau)} -- {line2(t0prime + tauprime + tprime_return)};

    # %\draw[->, thick, black!10!yellow] {line1(t0)} -- {line2(t0prime + tau)};
    # %\draw[->, thick, black!10!yellow] {line2(t0prime + tau)} -- {line1(t0 + tauprime + t_return)};
    # %\draw[->, thick, black!10!yellow] {line2(t0prime)} -- {line1(t0 + tauprime)};
    # %\draw[->, thick, black!10!yellow] {line1(t0 + tauprime)} -- {line2(t0prime + tau + tprime_return)};
    
    # %\draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(tauprime, t0)};
    # %\draw[->, thick, black!10!yellow] {light_r_to_l(0, t0 + tauprime)} -- {light_r_to_l(t_return, t0 + tauprime)};
    # %\draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tau, t0prime)};
    # %\draw[->, thick, black!10!yellow] {light_l_to_r(0, t0prime + tau)} -- {light_l_to_r(tprime_return, t0prime + tau)};
    
    # %\draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(tau, t0)};
    # %\draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime + tau)} -- {light_r_to_l(t_return, t0prime + tau)};
    # %\draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tauprime, t0prime)};
    # %\draw[->, thick, black!10!yellow] {light_l_to_r(0, t0 + tauprime)} -- {light_l_to_r(tprime_return, t0 + tauprime)};

    
    # % Tests
    # \draw {line1(t0)} node[circle, red, fill] {{}};
    # \draw {line1(t)} node[circle, red, fill] {{}};
    # \draw {line1(t_return)} node[circle, red, fill] {{}};
    # \draw {line2(t0prime)} node[circle, red, fill] {{}};
    # \draw {line2(tprime)} node[circle, red, fill] {{}};
    # \draw {line2(tprime_return)} node[circle, red, fill] {{}};




def moving_clocks_v2(v1, v2, tstart, tend, t0, t0prime=None):
    """
    Draw clock diagram for two moving observers which are at rest with
    respect to each other.

    If t1prime is None, we choose it such that clocks are synchronized.
    """
    v_avg = (v1 + v2) / 2
    # v_rel = np.abs(v1 - v2)

    def line1(t):  # Observer 1
        return (v1 * t, t)
        # return tuple([-1, 0] + [v1, 1] / np.sqrt(v1**2 + 1) * t)
    
    def line2(t):  # Observer 2
        return (v2 * t, t)
        # return tuple([1, 0] + [v2, 1] / np.sqrt(v2**2 + 1) * t)
    
    def line3(t):  # Referee
        # if v_avg != 0:
        #     return (v1 * v2 / v_avg * t, t)
        #     # return (t, t * v_avg / v1 / v2)
        # else:
        #     return (t, 0)

        # With formula from Mathematica, input was:
        # Solve[Norm[{v3*t, t}/Norm[{v3*t, t}] - {v2*t, t}/Norm[{v2*t, t}]] == 
        # Norm[{v3*t, t}/Norm[{v3*t, t}] - {v1*t, t}/Norm[{v1*t, t}]], v3, Reals, 
        # Assumptions ->  v1 \[Element] Reals && v2 \[Element] Reals && v3 \[Element] Reals && t \[Element] Reals]
        # if v_avg < 0:
        #     v3 = (-1 + v1 * v2) / (v1 + v2) - np.sqrt(1 + v1**2 + v2**2 + v1**2 * v2**2) / np.abs(v1 + v2)
        #     # print(v3)
        #     return (v3 * t, t)
        # elif v_avg > 0:
        #     v3 = (-1 + v1 * v2) / (v1 + v2) + np.sqrt(1 + v1**2 + v2**2 + v1**2 * v2**2) / np.abs(v1 + v2)
        #     # print(v3)
        #     return (v3 * t, t)
        # else:
        #     return (0, t)
        
        return (v_avg * t, t)

    def light_l_to_r(t, t0):
        return (v1 * t0 + t, t0 + t)

    def light_r_to_l(t, t0):
        return (v2 * t0 - t, t0 + t)
    

    # v1 /= np.sqrt(v1**2 + 1)
    # v2 /= np.sqrt(v2**2 + 1)
    

    if t0prime is None:
        # Set value such that clocks are synchronizeds
        t_ref = (v_avg * t0 - v1 * t0) / (1 - v_avg)  # Time where light from left observer reaches referee
        # tprime_ref = (1 + v2 * t0 + v1 * t0) / (1 + v2)  # Time where light from right observer reaches referee
        # tprime_ref = 1 / (1 + v2)  # Time where light from right observer reaches referee
        # t0prime = t0 + t_ref - tprime_ref
        t0prime = (t0 + t_ref - 1 / (1 + v_avg)) / (1 + (v2 - v_avg) / (1 + v_avg))
    

    # t = (2 + v2 * t0 - v1 * t0) / (1 - v2)  # t1 = t_-
    # # t = 2 / (1 - v2)  # t1 = t_-
    # # t_return = (2 + v1 * (t0 + t) - v2 * (t0 + t)) / (1 + v2)  # = t_+
    # t_return = (2 + v1 * (t0 + t) + v2 * (t0 + t)) / (1 + v2)  # = t_+; Not sure why +v2 * t0prime, but works
    # # t_return = 2 / (1 + v2)  # = t_+
    # tprime = (2 + v1 * t0prime + v2 * t0prime) / (1 + v2)  # Not sure why +v2 * t0prime, but works
    # # tprime = 2 / (1 + v2)
    # tprime_return = (2 + v2 * (t0prime + tprime) - v1 * (t0prime + tprime)) / (1 - v2)
    # # tprime_return = 2 / (1 - v2)


    # t = (2 + v2 * t0 - v1 * t0) / (1 - v2)  # t1 = t_-
    # tprime = (2 + v1 * t0 + v2 * t0prime) / (1 + v2)  # Not sure why +v2 * t0prime, but works
    # t_return = (2 + v1 * (t0 + t) + v2 * (t0prime + tprime)) / (1 + v2)  # = t_+; Not sure why +v2 * t0prime, but works
    # tprime_return = (2 + v2 * (t0prime + tprime) - v1 * (t0 + t)) / (1 - v2)


    # t = (2 + v_rel * t0) / (1 - v2)  # t1 = t_-
    # # t = 2 / (1 - v2)  # t1 = t_-
    # # t_return = (2 + v1 * (t0 + t) - v2 * (t0 + t)) / (1 + v2)  # = t_+
    # t_return = (2 + 2 * v_avg * (t0 + t)) / (1 + v2)  # = t_+; Not sure why +v2 * t0prime, but works
    # # t_return = 2 / (1 + v2)  # = t_+
    # tprime = (2 + 2 * v_avg * t0prime) / (1 + v2)  # Not sure why +v2 * t0prime, but works
    # # tprime = 2 / (1 + v2)
    # tprime_return = (2 + v_rel * (t0prime + tprime)) / (1 - v2)


    # This here seems to work
    # t0, t0prime = t0prime, t0  # Testing
    tau = (v2 * t0 - v1 * t0) / (1 - v2)  # t0 = t_-
    tauprime = (-v1 * t0prime + v2 * t0prime) / (1 + v1)
    # tau, tauprime = tauprime, tau  # Not needed if we switch t0 and t'0
    t_return = (-v1 * (t0 + tau) + v2 * (t0 + tau)) / (1 + v1)  # = t_+; Not sure why +v2 * t0prime, but works
    tprime_return = (v2 * (t0prime + tauprime) - v1 * (t0prime + tauprime)) / (1 - v2)  # = t'_+
    # t_return, tprime_return = tprime_return, t_return  # Not needed if we switch t0 and t'0


    t = (t0 + (t0 + tau + t_return)) / 2
    tprime = (t0prime + (t0prime + tauprime + tprime_return)) / 2


    print(tau, t0 * t_return)  # Not equal... Should they be?
    print(tauprime, t0prime * tprime_return)  # Not equal... Should they be?
    print(np.sqrt(1 - (v1 - v2)**2) * t, np.sqrt(1 - (v1 - v2)**2) * tprime)
    print(1 / np.sqrt(1 - (v1 - v2)**2) * t, 1 / np.sqrt(1 - (v1 - v2)**2) * tprime)

    
    return fr'''
\subfloat[Time $\tau'$ measured by observer $\mathcal{{O}}$]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(tau, t0)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0 + tau)} -- {light_r_to_l(t_return, t0 + tau)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line1(t0)} node[left] {{$t_-$}};
    \draw {line1(t0 + tau + t_return)} node[left] {{$t_+$}};
    \draw {line2(t0 + tau)} node[right] {{$\tau'$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
\end{{tikzpicture}}
}}\hspace*{{0.025\textwidth}}
%
\subfloat[Time $\tau$ measured by observer $\mathcal{{O}}'$]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tauprime, t0prime)};
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0prime + tauprime)} -- {light_l_to_r(tprime_return, t0prime + tauprime)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {line2(t0prime)} node[right] {{$t'_-$}};
    \draw {line2(t0prime + tauprime + tprime_return)} node[right] {{$t'_+$}};
    \draw {line1(t0prime + tauprime)} node[left] {{$\tau$}};
\end{{tikzpicture}}
}}\hspace*{{0.025\textwidth}}
%
\subfloat[Comparison of both observers]{{
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};
    \draw[->, thick] {line3(tstart)} -- {line3(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(tau, t0)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0 + tau)} -- {light_r_to_l(t_return, t0 + tau)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tauprime, t0prime)};
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0prime + tauprime)} -- {light_l_to_r(tprime_return, t0prime + tauprime)};
    
    %\draw[->, thick, black!10!yellow] {line1(t0)} -- {line2(tau)};
    %\draw[->, thick, black!10!yellow] {line2(tau)} -- {line1(t_return)};
    %\draw[->, thick, black!10!yellow] {line2(t0prime)} -- {line1(tauprime)};
    %\draw[->, thick, black!10!yellow] {line1(tauprime)} -- {line2(tprime_return)};

    %\draw[->, thick, black!10!yellow] {line1(t0)} -- {line2(t0prime + tauprime)};
    %\draw[->, thick, black!10!yellow] {line2(t0prime + tauprime)} -- {line1(t0 + tau + t_return)};
    %\draw[->, thick, black!10!yellow] {line2(t0prime)} -- {line1(t0 + tau)};
    %\draw[->, thick, black!10!yellow] {line1(t0 + tau)} -- {line2(t0prime + tauprime + tprime_return)};

    %\draw[->, thick, black!10!yellow] {line1(t0)} -- {line2(t0prime + tau)};
    %\draw[->, thick, black!10!yellow] {line2(t0prime + tau)} -- {line1(t0 + tauprime + t_return)};
    %\draw[->, thick, black!10!yellow] {line2(t0prime)} -- {line1(t0 + tauprime)};
    %\draw[->, thick, black!10!yellow] {line1(t0 + tauprime)} -- {line2(t0prime + tau + tprime_return)};
    
    %\draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(tauprime, t0)};
    %\draw[->, thick, black!10!yellow] {light_r_to_l(0, t0 + tauprime)} -- {light_r_to_l(t_return, t0 + tauprime)};
    %\draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tau, t0prime)};
    %\draw[->, thick, black!10!yellow] {light_l_to_r(0, t0prime + tau)} -- {light_l_to_r(tprime_return, t0prime + tau)};
    
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0)} -- {light_l_to_r(tau, t0)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime + tau)} -- {light_r_to_l(t_return, t0prime + tau)};
    \draw[->, thick, black!10!yellow] {light_r_to_l(0, t0prime)} -- {light_r_to_l(tauprime, t0prime)};
    \draw[->, thick, black!10!yellow] {light_l_to_r(0, t0 + tauprime)} -- {light_l_to_r(tprime_return, t0 + tauprime)};

    % Draw lines of simultaneity
    \draw[thick, blue] {line1(t)} -- {line2(tprime)};
    \draw[thick, blue] {line1(t0prime + tauprime)} -- {line2(t0 + tau)};

    % Make labels
    \draw {line1(tstart)} node[left] {{$\mathcal{{O}}$}};
    \draw {line2(t0 + tau)} node[right] {{$\tau'$}};
    \draw {line2(tprime)} node[right] {{$t'$}};
    \draw {line2(tstart)} node[right] {{$\mathcal{{O}}'$}};
    \draw {line1(t0prime + tauprime)} node[left] {{$\tau$}};
    \draw {line1(t)} node[left] {{$t$}};
    \draw {line3(tstart)} node[below] {{$\mathcal{{R}}$}};

    % Tests
    \draw {line1(t0)} node[circle, red, fill] {{}};
    \draw {line1(t)} node[circle, red, fill] {{}};
    \draw {line1(t_return)} node[circle, red, fill] {{}};
    \draw {line2(t0prime)} node[circle, red, fill] {{}};
    \draw {line2(tprime)} node[circle, red, fill] {{}};
    \draw {line2(tprime_return)} node[circle, red, fill] {{}};
\end{{tikzpicture}}
}}
    '''


def abb21(slope1, slope2, tstart, tend, t1, t2):
    def line1(t):
        return (slope1 * t, t)
    
    def line2(t):
        return (slope2 * t, t)
    
    def light(t0, t):
        return (slope1 * t0 + t, t0 + t)
        # \draw[->, thick, black!10!yellow] {light(t1, 0)} -- {light(t1, t1intersect)};
        # \draw[->, thick, black!10!yellow] {light(t2, 0)} -- {light(t2, t2intersect)};
    

    # Have to account for relative velocity, but also travel time of light (c=1)
    # slope1 * t1 + t_intersect = slope2 * t_intersect
    # <=> t_intersect = slope1 * t1 / (slope2 - 1)
    t1intersect = slope1 * t1 / (slope2 - 1)#t1 + (slope1 + slope2 + 1) * t1
    t2intersect = slope1 * t2 / (slope2 - 1)#t2 + (slope1 + slope2 + 1) * t2

    # slope1 * t1 + (t_intersect - t1) = slope2 * (t_intersect - t1)
    # <=> t_intersect = slope1 * t1 / (slope2 - 1)
    # t1intersect = (slope1 - 1 + slope2) * t1 / (slope2 - 1)#t1 + (slope1 + slope2 + 1) * t1
    # t2intersect = (slope1 - 1 + slope2) * t2 / (slope2 - 1)#t2 + (slope1 + slope2 + 1) * t2

    # t1intersect = (slope1 + slope2) * t1 / (slope2 - 1)#t1 + (slope1 + slope2 + 1) * t1
    # t2intersect = (slope1 + slope2) * t2 / (slope2 - 1)#t2 + (slope1 + slope2 + 1) * t2


    return fr'''
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {line1(t1)} -- {line2(t1intersect)};
    \draw[->, thick, black!10!yellow] {line1(t2)} -- {line2(t2intersect)};

    % Make labels
    \draw (0, {tstart}) node[below] {{$\mathcal{{O}}$}};
    \draw {line1(tend)} node[left] {{$\mathcal{{U}}$}};
    \draw {line1(t1)} node[left] {{$t_\mathcal{{U}}$}};
    \draw {line1(t2)} node[left] {{$t'_\mathcal{{U}}$}};
    \draw {line2(tend)} node[right] {{$\mathcal{{B}}$}};
    \draw {line2(t1intersect)} node[right] {{$t_\mathcal{{B}}$}};
    \draw {line2(t2intersect)} node[right] {{$t'_\mathcal{{B}}$}};
\end{{tikzpicture}}
    '''



def abb22(slope1, slope2, tstart, tend, t1):
    def line1(t):
        return (slope1 * t, t)
    
    def line2(t):
        return (slope2 * t, t)
    
    slope3 = (slope1 + slope2) / 2
    def line3(t, toffset=0):  # Referee
        return (slope3 * t, t + toffset)
    
    def light_positive(t0, t):
        return (slope1 * t0 + t, t0 + t)
    
    def light_negative(t0, t):
        return (slope2 * t0 - t, t0 + t)
    

    # t1intersect = np.sqrt(1 - (slope1)**2) * t1
    # t2intersect = np.sqrt(1 - (slope2)**2) * t1
    tau_u = np.sqrt(1 - (slope1 + slope2)**2) * t1
    tau_b = t1

    # slope1 * t1 + t_intersect = slope3 * t_intersect
    # <=> t_intersect = slope1 * t1 / (slope3 - 1)
    t_intersect_1 = slope1 * tau_u / (slope3 - 1)
    # slope2 * t1 - t_intersect = slope3 * t_intersect
    # <=> t_intersect = slope2 * t1 / (slope3 + 1)
    t_intersect_2 = slope2 * tau_b / (slope3 - 1)
    t_intersect_2 = t_intersect_1
    # slope3 * tintersect - t_intersect = slope2 * t1
    # <=> t_intersect = slope2 * t1 / (slope3 - 1)
    # t_intersect_2 = slope2 * tau_b / (slope3 - 1)

    # slope1 * t1 + (t_intersect - t1) = slope3 * (t_intersect - t1)
    # <=> t_intersect = (slope1 - 1 + slope3) * t1 / (slope3 - 1)
    # t_intersect_1 = (slope1 - 1 + slope3) * tau_u / (slope3 - 1)#t1 + (slope1 + slope3 + 1) * t1
    # slope1 * t1 - t_intersect = slope3 * t_intersect
    # <=> t_intersect = slope1 * t1 / (slope3 + 1)
    # t_intersect_2 = slope1 * tau_u / (slope3 + 1)#t1 + (slope1 + slope3 - 1) * t1
    # slope1 * t1 - (t_intersect - t1) = slope3 * (t_intersect - t1)
    # <=> t_intersect = slope1 * t1 / (slope3 + 1)
    # t_intersect_2 = (slope1 + 1 + slope3) * tau_u / (slope3 + 1)#t1 + (slope1 + slope3 + 1) * t1


    # # slope1 * tau_u + t_intersect = slope2 * tau_b - t_intersect
    # t_intersect_1 = (- slope1 * tau_u + slope2 * tau_b) / 2
    # t_intersect_2 = (slope1 * tau_u - slope2 * tau_b) / 2
    


    return fr'''
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers and referee
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};
    \draw[->, thick] {line3(tstart)} -- {line3(tend)};

    % Draw line of simultaneity
    \draw[thick, black!10!green] {line1(tau_u)} -- {line2(tau_b)};

    % Draw parallelogram
    \draw[thick] {line1(tau_u)} -- {line3(t_intersect_1, toffset=tau_u)};
    \draw[thick] {line1(tau_u)} -- {line3(t_intersect_2, toffset=tau_u - 2 * t_intersect_2)};
    \draw[thick] {line3(t_intersect_1, toffset=tau_u)} -- {line2(tau_b)};
    \draw[thick] {line3(t_intersect_2, toffset=tau_u - 2 * t_intersect_2)} -- {line2(tau_b)};
    %\draw[thick] {light_negative(tau_b, 0)} -- {light_negative(tau_b, t_intersect_1)};
    %\draw[thick] {light_positive(tau_u, 0)} -- {light_positive(tau_u, t_intersect_2)};

    % Make labels
    \draw (0, {tstart}) node[below] {{$\mathcal{{O}}$}};
    \draw {line1(tend)} node[left] {{$\mathcal{{U}}$}};
    \draw {line1(tau_u)} node[left] {{$\tau_\mathcal{{U}}$}};
    \draw {line3(tend)} node[left] {{$\mathcal{{S}}$}};
    \draw {line2(tend)} node[right] {{$\mathcal{{B}}$}};
    \draw {line2(tau_b)} node[right] {{$\tau_\mathcal{{B}}$}};
\end{{tikzpicture}}
    '''



def abb23(slope1, slope2, tstart, tend, t1, t2):
    def line1(t):
        return (slope1 * t, t)
    
    def line2(t):
        return (slope2 * t, t)
    
    slope3 = (slope1 + slope2) / 2
    def line3(t):  # Referee
        return (slope3 * t, t)
    

    t1intersect = t1 + (slope1 + slope2 + 1) * t1
    t2intersect = t2 + (slope1 + slope2 + 1) * t2


    return fr'''
\begin{{tikzpicture}}[thick, >={{[inset=0,angle'=27]Stealth}}]
    % Draw lines of observers
    \draw[->, thick] {line1(tstart)} -- {line1(tend)};
    \draw[->, thick] {line2(tstart)} -- {line2(tend)};

    % Draw trajectories of light
    \draw[->, thick, black!10!yellow] {line1(t1)} -- {line2(t1intersect)};
    \draw[->, thick, black!10!yellow] {line1(t2)} -- {line2(t2intersect)};

    % Make labels
    \draw (0, {tstart}) node[below] {{$\mathcal{{O}}$}};
    \draw {line1(tend)} node[left] {{$\mathcal{{U}}$}};
    \draw {line1(t1)} node[left] {{$t_\mathcal{{U}}$}};
    \draw {line1(t2)} node[left] {{$t'_\mathcal{{U}}$}};
    \draw {line2(tend)} node[right] {{$\mathcal{{B}}$}};
    \draw {line2(t1intersect)} node[right] {{$t_\mathcal{{B}}$}};
    \draw {line2(t2intersect)} node[right] {{$t'_\mathcal{{B}}$}};
\end{{tikzpicture}}
    '''



if __name__ == '__main__':
    # print(resting_clocks(0, 5, 0.5))  # Good

    # print(moving_clocks_resting_wrt_each_other(0.3, 0, 6, 0.5))  # Good

    # For comparison of simultaneous etc
    print(resting_clocks_single_fig(0, 6, 0.5))
    print(moving_clocks_resting_wrt_each_other_single_fig(0.3, 0, 6, 0.5))
    # print(resting_clocks(0, 6, 0.5, x_shift=0.3 * (6 - 0.5) / 2 * 1 / np.sqrt(1 - 0.3**2), x_distance=1 / np.sqrt(1 - 0.3**2)))
    print(resting_clocks_fitted(0, 6, 0.5, 0.3))

    # print(moving_clocks(0.1, 0.4, 0, 6, 0.5))
    # print(moving_clocks(-0.1, 0.4, -1, 7, -0.5))  # Does not work with negative velocities yet -> does now
    # print(moving_clocks(-0.3, 0.3, -1.5, 5, -1))  # Good
    # print(moving_clocks(-0.2, 0.4, -4, 4, -2, -2 * np.sqrt(1 - 0.4**2)))
    # Multiplying with that looks best, but don't know why because it is wrong velocity...
    # Dividing by np.sqrt(1 - 0.4**2) also does not help; in fact, -2 for both looks best (but not better than automatic)
    # print(moving_clocks(-0.1, 0.4, -4, 3, -2))  # This makes most sense, only then lines between tau, tau' and t, t' parallel (needed!)
    # print(moving_clocks(-0.1, 0.4, 0.0, 7, 2, 2 / np.sqrt(1 - 0.5**2) - 0.2))
    # print(moving_clocks(0.0, 0.5, 0.0, 7, 2, 2 / np.sqrt(1 - 0.5**2) - 0.1))
    # print(moving_clocks(0.0, 0.5, -4, 4, -2))
    # print(moving_clocks(0.1, 0.4, -1, 7, 0.5, np.sqrt(1 - 0.3**2) * 0.5))
    # print(moving_clocks(0.1, 0.4, -1, 7, 0.5))

    # print(abb21(-0.4, 0.75, -0.5, 5, 1, 2))

    # print(abb22(-0.4, 0.75, -0.5, 5, 2))

    # with open('../pictures/abb21.tex', 'w') as file:
    #     file.write(abb21(-0.5, 0.5, -0.5, 5, 1, 2))