## Introduction

Simulation code for proof of concept adaptive spread-spectrum underwater communications using auxiliary-vector filtering. Includes an [underwater channel simulator](http://millitsa.coe.neu.edu/?q=projects) developed by P. Qarabaqi and M. Stojanovic at Northeastern to generate time-varying channel realizations and a [library to generate auxiliary-vector sequences](http://www.eng.buffalo.edu/~pmarkopo/av.php) from P. Markopoulos at SUNY Buffalo. 

This software is free to use and modify (see LICENSE). If helpful in your own work, please cite 

```
@inproceedings{gannon2018,
	title={Short Data Record Filtering for Adaptive Underwater Acoustic Communications},
	author={Gannon, Adam and Balakrishnan, Sarankumar and Sklivanitis, George and Pados, Dimitris A. and Batalama, Stella N.},
	booktitle={IEEE Sensor Array and Multichannel Signal Processing Workshop (SAM)},
	year={2018},
}
```

## Dependencies
Requires MATLAB. 

## Scripts 

### Main Scripts

#### full_loop_ct.m
Simulates the adaptive code length feedback system. Packets propagate over a channel generated from the acoustic_channel_simulator package. Code length is determined based on instantaneous channel gain and fed back to the transmitter with a simulated propagation delay. 

#### postprocessing_full.m
Plots additional comparisons from full_loop_ct.m

### Additional Scripts

#### ber_vs_sample_support.m
Calculates the bit error rate performance of each filter (MF, MVDR, AV) as a function of sample support. 

#### ber_vs_sinr_av.m
For each code length, calculates the average bit error rate as a function of transmit signal energy. 

#### sinr_vs_chan_gain.m
Calculates the performance of each code length as a function of the channel gain. 


