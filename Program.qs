namespace hamiltonian_simulation {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Oracles;
    open auxiliary;
    

    @EntryPoint()
    operation SayHello() : Unit {
        Message("Hello quantum world!");
    }

    // Welcome to the introduction to simulation coding assignment!
    // Written by: Christopher Kang
    //
    // Learning goals
    // -> Understand the core fundamentals of simulation
    // -> Understand which gates are necessary for simulation
    // -> Implement and understand the Variational Quantum Eigensolver
    // -> Implement and understand Trotterization

    // ---------------------------------------------------
    // I. FUNDAMENTALS OF SIMULATION
    // ---------------------------------------------------
    // 
    // Suppose we have a Hamiltonian H that we'd like to emulate.
    // Namely, we'd like to implement $ e^{-iHt} $ as a unitary operation
    // upon some starting quantum state $ |\psi> $, representing a starting
    // molecular state. How would we implement $ e^{-iHt} $?
    //
    // First, we can start with a simple Hamiltonian: H = Z
    // So, we'd like to implement $ e^{-iHt} = e^{-iZt} $.
    // Let's go ahead and do that (reference the Q# documentation for the gate)
    operation SimpleHamiltonian1(q : Qubit) : Unit {
        Exp([PauliZ], 1.0, [q]);
    }

    // We could also try a more complicated Hamiltonian, one which acts on multiple qubits
    // Let's try H = X \kron Z (X on qubit 0, Z on qubit 1).
    operation SimpleHamiltonian2(q : Qubit[]) : Unit {
        Exp([PauliX], 1.0, [q[0]]);
        Exp([PauliZ], 1.0, [q[1]]);
    }

    // We'd also need to prepare suitable states for simulation
    // Let's say that our starting molecule is \sqrt{2/3}|01> + \sqrt{1/3}|10>
    // Can you prepare the state? (Assume that the ancilla is given in 00)
    // (PS: this state could represent an electron in a lower orbital with 2/3 probability
    // and simultaneously in the upper orbital with 1/3 probability)
    operation PrepareState(q : Qubit[]) : Unit {
        // 
    }

    // WRITTEN QUESTION 1:
    // WHAT IS THE MATRIX REPRESENTATION OF THE FOLLOWING HAMILTONIAN?
    // Z \KRON X + X \KRON Z

    // ---------------------------------------------------
    // II. Exploration of Trotterization
    // ---------------------------------------------------
    // 
    // As seen above, not all Hamiltonians can be easily implemented.
    // Worse, we're not even guaranteed the Hamiltonian is unitary,
    // only that it is hermitian (which leads to the exponentiated matrix
    // being unitary). How can we get around this?
    //
    // We leverage the approach known as Trotterization. Further information
    // is provided in the handout, but the general idea is that
    // we may split the Hamiltonian into easily simulatable units
    // with minimal error. For example, taking the example above:
    // e^{-i(Z \kron X + X \kron Z)t} \approx e^{-i(Z \kron X)t} e^{-i(X\kronZ)t} + O(t^2m)
    //
    // Let's try to implement this directly.
    operation HamiltonianLeft(q : Qubit[]) : Unit is Adj + Ctl {
        // IMPLEMENT ME! (❁´◡`❁)
        // (Z \kron X)
        Exp([PauliZ, PauliX], 1.0, q);
    }

    operation HamiltonianRight(q : Qubit[]) : Unit is Adj + Ctl{
        // IMPLEMENT ME! (❁´◡`❁)
        Exp([PauliX, PauliZ], 1.0, q);
    }

    operation TotalHamiltonian(q : Qubit[]) : Unit is Adj + Ctl{
        HamiltonianLeft(q);
        HamiltonianRight(q);
    }

    operation StatePrepRoutine(q : Qubit[]) : Unit {
        // For simplicity, let's use the equal superposition over |01>, |10>
        // IMPLEMENT ME!
        H(q[0]);
        CNOT(q[0], q[1]);
        X(q[0]);
    }

    operation PrimitiveTrotterization() : Double {
        // YOU DO NOT NEED TO MODIFY THIS CODE!
        // IT IS A DEMO ON HOW TO GET THE PHASE
        mutable phaseEst = 0.0;
        let PRECISION = 7;

        using ((arr, anc) = (Qubit[2], Qubit[PRECISION])) {
            let precision = BigEndian(anc);
            StatePrepRoutine(arr);
            QuantumPhaseEstimation(OracleToDiscrete(TotalHamiltonian), arr, precision);
            let result = MeasureInteger(BigEndianAsLittleEndian(precision));
            set phaseEst = IntAsDouble(result) / IntAsDouble(1 <<< PRECISION);
        }

        return phaseEst;
    }

    // While this estimate may have been okay, we can do better.
    // Refer to the guide on the exact scaling of Trotterization
    // We'll now implement Trotterization with trotter step 0.5
    // So, we'll need to implement the half steps each twice
    // Or, (e^{-i (Z \kron X) t / 2} e^{-i (X \kron Z) t / 2})^2
    operation HalfHamiltonianLeft(q : Qubit[]) : Unit is Adj+Ctl {
        //
        Exp([PauliZ, PauliX], 0.5, q);
    }

    operation HalfHamiltonianRight(q : Qubit[]) : Unit is Adj+Ctl {
        //
        Exp([PauliX, PauliZ], 0.5, q);
    }

    operation HalfTotalHamiltonian(q : Qubit[]) : Unit is Adj+Ctl {
        for (_ in 0..2) {
            HalfHamiltonianLeft(q);
            HalfHamiltonianRight(q);
        }
    }

    operation DoubleStepTrotterization() : Double {
        return EstimatePhase(HalfTotalHamiltonian, StatePrepRoutine, 7);
    }

    // Finally, go ahead and implement a general purpose version of Trotterization
    // Use the provided signature. The number of steps used is 
    operation FullTrotterization(oracles : ((Qubit[], Double) => Unit is Adj+Ctl)[], steps : Int) : Double {
        //
        Double stepSize = 1.0 / IntAsDouble(steps); 

        using (arr = Qubit[2]) {
            //
        }

        return EstimatePhase(_FullTrotterOracle(_, oracles, 5), NoOp<Qubit[]>(_), 7);
    }

    operation _FullTrotterOracle(reg : Qubit[], oracles : ((Qubit[], Double) => Unit is Adj+Ctl)[], steps : Int) : Unit is Adj+Ctl {
        //
        for (_ in 0..steps) {
            let stepSize = 1.0 / IntAsDouble(steps);
            for (oracle in oracles) { 
                oracle(reg, stepSize);
            }
        }
    }

    operation _TrotterStep(reg : Qubit[], hamiltonian : (Qubit[] => Unit), stepSize : Double) : Unit {
        //
    }


    // ---------------------------------------------------
    // III. Exploration of VQEs
    // ---------------------------------------------------
    //
    // Trotterization is an incredibly powerful technique because its error bound
    // scales well with the Hamiltonian size and desired timestep. Unfortunately,
    // Trotterization will still be difficult to implement on near-term devices because
    // it often requires circuits that are hundreds, if not thousands, of operations deep.
    // Moreover, it requires Quantum Phase Estimation, which can consume valuable qubits.
    //
    // The Variational Quantum Eigensolver (VQE) was developed as a near-term alternative.
    // The "Variational" part comes from the fact that the state preparation routine
    // can be parametrized and then optimized classically. The algorithm is an "Eigensolver"
    // because finding the minimum energy level is essentially eigenvalue finding.
    //
    // The core insight with the VQE is that the Hamiltonian computation doesn't need to all
    // happen in the quantum computer. Consider the following Hamiltonian:
    //
    // H = 0.01 * X \kron Z 
    //
    // In a typical quantum computer, implementing this Hamiltonian would be incredibly difficult;
    // We lack the precision to implement such small angles for these components
    // However, we can easily find the phase induced by X \kron Z, and simply exponentiate the values
    // and multiply them by their coefficients.
    //
    // This property is linearity; namely:
    // IF H = H_1 + H_2 + ... + H_n and H|\psi> = E|\psi>, then E can be approximated as H_1|\psi> + ..

    operation ImplementVQE(hamiltonian : (Qubit[] => Unit)[], size : Int) : Double {
        Double final = 0.0;
        Int NUMBER_OF_ITERATIONS = 10;
        phase_est_alg = QuantumPhaseEstimation();

        using (arr = Qubit[size]) {
            for (sub : hamiltonian) {
                Double total = 0.0;
                for (iter = 0; iter < NUMBER_OF_ITERATIONS; iter++) {
                    StatePrep(arr);
                    total += phase_est_alg(sub);
                    ResetAll(arr);
                }
                total = total / NUMBER_OF_ITERATIONS;
                final += total;
            }
        }

        return final;
    }



    // ---------------------------------------------------
    // IV. It's your favorite molecule...WATER!
    // ---------------------------------------------------
    //
    // The code you've made is awesome, and we're going to take it a step further.
    // The main use of chemical simulations is to identify energy amounts and
    // molecular structure - if we've found the "ground" or stable state, the 
    // molecular energy will be minimized.
    //
    // We are using data from the NWChem package, a computational chemistry
    // package developed by PNNL. This produces some data for us to use.
    // We will abstract a lot of this away; if you are interested in how it works,
    // Microsoft has excellent documents on Second Quantization and the 
    // Jordan-Wigner transformation [INSERT LINK HERE!]
    //
    // If your code works, a few things should happen:
    // 1. You will be able to identify the most "stable" of the molecules we provide
    // (and tell us which one is most like the water we drink!)
    // 2. You will be able to compare accuracy between VQE and Trotterization, especially
    // as the number of steps for Trotterization increases + as the number of queries for 
    // VQE increases.
    // 3. [SOMETHING ELSE]

    WATER_MOLECULE_1 = 0;
    WATER_MOLECULE_2 = 0;
    WATER_MOLECULE_3 = 0;

    operation EstimateWaterEnergy() : Unit {
        //
    }

    operation CompareTrotterizationToVQE() : Unit {
        //
    }

    // ---------------------------------------------------
    // V. [Extra Credit] Implementing gates for Trotterization
    // ---------------------------------------------------
    //
    // It's easy for us to abstract away the gates necessary for Trotterization;
    // we simply say "Exp" to exponentiate the gate. However, what's the actual 
    // implementation? Especially if we are limited to a typical gateset.
    //
    // Question 1: What is the matrix of e^{-iZ}?
    //
    // Question 2: How does e^{-iZ} act on |0>? On |1>?
}
