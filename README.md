## Introduction

Simulation code for proof of concept adaptive spread-spectrum underwater communications using auxiliary-vector filtering. Includes an [underwater channel simulator](http://millitsa.coe.neu.edu/?q=projects) developed by P. Qarabaqi and M. Stojanovic at Northeastern to generate time-varying channel realizations and a [library to generate auxiliary-vector sequences](http://www.eng.buffalo.edu/~pmarkopo/av.php) from P. Markopoulos at SUNY Buffalo.

This software is free to use and modify (see LICENSE). Please include credit to the above authors of the included software libraries and maintain links to their websites in any derivative work.

If helpful in your own research, please cite
```
@inproceedings{gannon2018,
	title={Short Data Record Filtering for Adaptive Underwater Acoustic Communications},
	author={Gannon, Adam and Balakrishnan, Sarankumar and Sklivanitis, George and Pados, Dimitris A. and Batalama, Stella N.},
	booktitle={IEEE Sensor Array and Multichannel Signal Processing Workshop (SAM)},
	month={july},
	address={Sheffield, England}
	year={2018},
}
```


## Dependencies
Requires MATLAB. Tested on 2011b and 2017a. The simulation will use the DSP and Communications toolboxes, if available, but should also work without it.

## Scripts

### Preliminary Scripts

#### acoustic_channel_simulator/channel_simulator.m
First, run the Qarabaqi/Stojanovic channel simulator to generate the channel from the parameter files provided in ``acoustic_channel_simulator/channel_simulator/channels/``.

### Main Scripts

#### full_loop_ct.m
Simulates the adaptive code length feedback system. Packets propagate over a channel generated from the acoustic_channel_simulator package. Code length is determined based on instantaneous channel gain and fed back to the transmitter with a simulated propagation delay.

#### postprocessing_full.m
Plots additional comparisons from ``full_loop_ct``. The script assumes that ``full_loop_ct`` has been run once without adaptation and once with adaptation to generate two .mat files which are read in.

### Additional Scripts

#### ber_vs_sample_support.m
Calculates the bit error rate performance of each filter (MF, MVDR, AV) as a function of sample support.

#### ber_vs_sinr_av.m
For each code length, calculates the average bit error rate as a function of transmit signal energy.

#### sinr_vs_chan_gain.m
Calculates the performance of each code length as a function of the channel gain.
