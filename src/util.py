import math
import matplotlib.pyplot as plt
import numpy as np
from typing import *

from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram, circuit_drawer
import qiskit.visualization
from qiskit.quantum_info import Statevector, Operator, DensityMatrix


# -----------------------------------------------------------------------------
# Binary
# -----------------------------------------------------------------------------

def to_binary(num_qubits: int, j: int) -> list[int]:
    return [int(x) for x in f"{j:0{num_qubits}b}"]


def idx_to_basis(num_qubits: int, j: int) -> Statevector:
    vecs = [zero if x == 0 else one for x in to_binary(num_qubits, j)]
    acc = vecs[-1]
    for v in vecs[-2::-1]:
        acc = v ^ acc
    return acc


# -----------------------------------------------------------------------------
# Measurement
# -----------------------------------------------------------------------------

def Pretty(arr: np.ndarray):
    return Operator(arr).draw("latex")


# -----------------------------------------------------------------------------
# Gates
# -----------------------------------------------------------------------------

def mk_X():
    qc = QuantumCircuit(1)
    qc.x(0)
    return Operator(qc)


def mk_Y():
    qc = QuantumCircuit(1)
    qc.y(0)
    return Operator(qc)


def mk_Z():
    qc = QuantumCircuit(1)
    qc.z(0)
    return Operator(qc)


def mk_H():
    qc = QuantumCircuit(1)
    qc.h(0)
    return Operator(qc)


def mk_CNOT():
    qc = QuantumCircuit(2)
    qc.cx(0, 1)  # 0 is control, 1 is target
    return Operator(qc)

X = mk_X()
Y = mk_Y()
Z = mk_Z()
H = mk_H()
CNOT = mk_CNOT()


# -----------------------------------------------------------------------------
# Qubits
# -----------------------------------------------------------------------------

zero = Statevector(np.array([1.0 + 0j, 0j]))
one = Statevector(np.array([0j, 1.0 + 0j]))


def qubit_condition(q: Union[np.array, Statevector]) -> np.array:
    return q[0]*np.conjugate(q[0]) + q[1]*np.conjugate(q[1])


def is_qubit(q: Union[np.array, Statevector]) -> bool:
    return np.allclose(np.array([1.]), qubit_condition(q))


def demonstrate_measure(q):
    sim = AerSimulator()
    
    # Don't worry about this code for now
    qc = QuantumCircuit(1, 1)
    qc.initialize(q, 0)        
    qc.measure(0, 0)
    
    results = sim.run(qc, shots=10).result()
    answer = results.get_counts()
    return plot_histogram(answer)


# -----------------------------------------------------------------------------
# Bloch Sphere
# -----------------------------------------------------------------------------

def statevector_to_bloch_vector(state):
    bloch_x = 2 * np.real(state[0] * np.conj(state[1]))
    bloch_y = 2 * np.imag(state[1] * np.conj(state[0]))
    bloch_z = np.abs(state[0])**2 - np.abs(state[1])**2
    return [bloch_x, bloch_y, bloch_z]


def plot_bloch_vector(q: Union[np.ndarray, Statevector], ax=None, title=None) -> None:
    return qiskit.visualization.plot_bloch_vector(statevector_to_bloch_vector(q), ax=ax, title=title)


class PlotGateOpOnBloch:
    def __enter__(self):
        fig = plt.figure(figsize=(8, 4))
        self.ax1 = fig.add_subplot(121, projection='3d')
        self.ax2 = fig.add_subplot(122, projection='3d')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass


# -----------------------------------------------------------------------------
# Density Matrix
# -----------------------------------------------------------------------------

def measure_outcome(rho: DensityMatrix, M_m: np.ndarray):
    # Probability of measuring m
    p_m = np.trace(np.conjugate(M_m.transpose()) @ M_m @ rho.data)
    # Resulting state after m
    rho_p = (M_m @ rho.data @ np.conjugate(M_m.transpose())) / p_m
    return p_m, rho_p


# -----------------------------------------------------------------------------
# Bits
# -----------------------------------------------------------------------------

def enum_bits(n: int) -> list[list[bool]]:
    if n == 1:
        return [[False], [True]]
    else:
        acc = []
        for x in enum_bits(n-1):
            acc += [[False] + x]
            acc += [[True] + x]
        return acc


# -----------------------------------------------------------------------------
# Number theory
# -----------------------------------------------------------------------------

def prime_seive(n):
    def seive(primes, i):
        for prime in primes:
            if prime > math.ceil(math.sqrt(n)):
                break
            if i % prime == 0:
                return False
        return True
    
    primes = [2, 3, 5]
    for i in range(6, n, 6):
        if seive(primes, i - 1):
            primes += [i - 1]
        if seive(primes, i + 1):
            primes += [i + 1]
    return primes


def euler_totient(N):
    cnt = 0
    for i in range(N):
        if math.gcd(i, N) == 1:
           cnt += 1
    return cnt


def factor2_to_order(N: int, a: int, s: int) -> int:
    # Check that s is even
    if s % 2 != 0:
        print(s)
        return None
    
    guesses = [math.gcd(a**(s//2)-1, N), math.gcd(a**(s//2)+1, N)]
    for guess in guesses:
        # Check to see if guess is a factor
        if guess not in [1, N] and (N % guess) == 0:
            return guess


# -----------------------------------------------------------------------------
# QFT
# -----------------------------------------------------------------------------

def qft(circuit: QuantumCircuit, n: int) -> QuantumCircuit:
    def go(circuit, n):
        if n == 0:
            return circuit
        else:
            n -= 1
            # Apply H
            circuit.h(n)
            # Apply CR
            for qubit in range(n):
                circuit.cp(np.pi/2**(n-qubit), qubit, n)
            # Recurrence relation
            return go(circuit, n)
    circuit = go(circuit, n) 

    # Take care of little endian
    for qubit in range(n//2):
        circuit.swap(qubit, n-qubit-1)
    return circuit


def qft_dagger(qc, n):
    # Swap for little endian
    for qubit in range(n//2):
        qc.swap(qubit, n-qubit-1)
    for j in range(n):
        # Inverse of controlled rotation is inverse rotation
        for m in range(j):
            qc.cp(-math.pi/float(2**(j-m)), m, j)
        # Hadamard is own inverse
        qc.h(j)


# -----------------------------------------------------------------------------
# Modular Exponentiation
# -----------------------------------------------------------------------------

def amod15(a, power):
    # Adapted From: https://qiskit.org/textbook/ch-algorithms/shor.html
    if a not in [2,4,7,8,11,13]:
        raise ValueError("'a' must be 2,4,7,8,11 or 13")
    U = QuantumCircuit(4)        
    for iteration in range(power):
        if a in [2,13]:
            U.swap(2,3); U.swap(1,2); U.swap(0,1)
        if a in [7,8]:
            U.swap(0,1); U.swap(1,2); U.swap(2,3)
        if a in [4, 11]:
            U.swap(1,3); U.swap(0,2)
        if a in [7,11,13]:
            for q in range(4):
                U.x(q)
    U = U.to_gate(); U.name = "%i^%i mod 15" % (a, power); 
    return U


# -----------------------------------------------------------------------------
# Quantum Phase Estimation
# -----------------------------------------------------------------------------

def qpe(f_qc_U, n_count):
    # Create and set up circuit
    reg1_size = f_qc_U(0).num_qubits
    qc = QuantumCircuit(n_count+reg1_size, n_count)
    # Step 1: Prepare our eigenstate |1>:
    qc.x(n_count)
    # Step 2: Apply H-Gates to counting qubits.
    qc.barrier()
    for qubit in range(n_count):
        qc.h(qubit)
    # Step 3: Do the controlled-U operations:
    for q in range(n_count):
        qc_cU = f_qc_U(2**q).control() # amod15(a, 2**q).control()
        qc.append(qc_cU, [q] + [n_count+i for i in range(reg1_size)])
    # Step 4: Perform inverse QFT
    qft_dagger(qc, 4)
    # Step 5: Measure first register
    qc.barrier()
    for n in range(n_count):
        qc.measure(n, n)
    return qc

# qc_qpe = qpe(lambda q: amod15(7, q), 4)


# -----------------------------------------------------------------------------
# Drawing vectors
# -----------------------------------------------------------------------------

def draw_vecs(vecs):
    V = []; origin_x = []; origin_y = []
    for x in vecs:
        if type(x) is tuple:
            V += [x[1]]; origin_x += [x[0][0]]; origin_y += [x[0][1]]
        else:
            V += [x]; origin_x += [0]; origin_y += [0]
    V = np.array(V)
    origin = np.array([origin_x, origin_y])
    W = V + np.transpose(origin)    
    fig = plt.figure(); ax = fig.add_subplot()
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.quiver(*origin, V[:,0], V[:,1], angles='xy', scale_units='xy', color=['r', 'g', 'b'], scale=1)
    x_scale = .1 * (np.max(W[:,0]) - min(0, np.min(W[:,0])))
    y_scale = .1 * (np.max(W[:,1]) - min(0, np.min(W[:,1])))
    l = min(min(0, np.min(W[:,0])) - x_scale, min(0, np.min(W[:,1])) - y_scale)
    r = max(np.max(W[:,0]) + x_scale, np.max(W[:,1]) + y_scale)
    plt.xlim(l, r); plt.ylim(l, r)