import numpy as np


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



print(abb21(-0.4, 0.75, -0.5, 5, 1, 2))

print(abb22(-0.4, 0.75, -0.5, 5, 2))

# with open('../pictures/abb21.tex', 'w') as file:
#     file.write(abb21(-0.5, 0.5, -0.5, 5, 1, 2))