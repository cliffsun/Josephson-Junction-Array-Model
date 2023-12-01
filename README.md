# Simulation-of-non-trivial-Josephson-Junction-SQUID

This code describes the behavior of the Critical Current in one SQUID (Superconducting Quantum Interference Devices) when subject to a quantized magnetic flux.

The array "arrayOfJunctions" signifies where each Josephson Junction is located (and thus where current is flowing). 

If the array is as follows: [0, 0.25, 0.5, 1] -> that means that current is flowing from 0 - 0.25 and 0.5 - 1. The space in between 0.25 and 0.5 shows where current _isn't_ flowing. 

This array must have an even number of values. The numbers in the array represent the percentage that the junction is taking up relative to the circuit. So from 0 - 0.25 means that this junction would take up 25% of the circuit (which is a lot!). 

Any questions about the code: you can contact me @cliffxs2@illinois.edu

Credit: Harshvardhan Mantry and Professor Alexey Bezryadin
