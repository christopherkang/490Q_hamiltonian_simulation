namespace auxiliary {
    open Microsoft.Quantum.Arithmetic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Characterization;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Oracles;


    operation EstimatePhase(oracle : (Qubit[] => Unit is Adj+Ctl), 
                            statePrepRoutine : (Qubit[] => Unit), precision : Int) : Double {
        mutable phaseEst = 0.0;
        using ((arr, anc) = (Qubit[2], Qubit[precision])) {
            // setup the groundwork
            let discrete = OracleToDiscrete(oracle);
            let bit_anc = BigEndian(anc);

            // apply the operations
            statePrepRoutine(arr);
            QuantumPhaseEstimation(discrete, arr, bit_anc);

            // prepare the result
            let result = MeasureInteger(BigEndianAsLittleEndian(bit_anc));
            set phaseEst = IntAsDouble(result) / IntAsDouble(1 <<< precision);
        }

        return phaseEst;
    }

    operation ParametrizedTrotterGate(reg : Qubit[], step : Double) : Unit is Adj+Ctl {
        //
    }

    operation SimpleHamiltonian(q : Qubit) : Unit {
        // This is a simple Hamiltonian for your testing purposes
        // It implements the following transformation:
        // H = Z
    }
}


