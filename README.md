# SyntheticECG
An algorithm for Generation of ECG Signals using Fourier Series Approximation.

This algorithm was presented in a paper titled "Generation of ECG for Heart Block Cases" in the Springer Journal.

Link to the paper - https://bit.ly/2QPEsuu

## Introduction

An electrocardiogram, short for ECG or EKG, is a tool used to measure the electrical
activities of the heart. This is measured by attaching a pair of electrodes on the
surface of the chest and on various other body parts, for accuracy. Small differences
in the potential as a consequence of cardiac muscle depolarization (contraction) and
immediate re-polarization (relaxation) constitutes for one complete heartbeat. Every
such heartbeat is measured and plotted on the graph.
Real-time ECG signals are contaminated with various noises and artefacts,
resulting in a signal which would pose difficulty in analyzing the output. Hence,
there is a requirement for a signal to be processed, in order to interpret an accurate
result from the graph which enables to infer required diagnostic information. Albeit,
to achieve such outcomes, a reliable signal processing is required, which consumes
time and ample efforts to come to a viable solution. Therefore, there is a need for
standardized reference. Additionally, accuracy is of not a major concern, if it comes
down to analyzing the data for academic studies.

This paper throws a light on generating artificial ECG signals using
geometrical features and Fourier series approximation, with the goal to simulate
several heart blocks. Furthermore, this study is one half of the entire project
undertaken to generate and detect various abnormal heartbeats.




## Note!

Algorithm.m is the basic foundation on which the enitre project is build. It gives an in sight to the algorithm for gnerating the ECG signals.
Normal Sinus Rythm (NSR) is generated as an example.
 
Block folder contains couple of heart block cases. 


