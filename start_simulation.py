#command [percent of slacks] [heterdesign, e.g. _GeSi without _
import os
import math
import string
import struct
import sys 

base_dir = './'
os.chdir(base_dir)
trials = 419430400
n_blocks = 128
n_treads_p_block = 128
initial_state = 1 # 1: ap2p, 0: p2ap

STTorPS = 0 # 0: STT, 1:Precessional Swiching

sense_time = 0
pulse_start = 1.5
pulse_end = 17
pulse_step = 0.5
exist_results = list()
finished_results = './finished_results.txt'
template_pulse_filename = './template_pulse.txt'
pulse_filename = './var_input.txt'
temperature = 300;
V_p = -0.7
V_ap = -0.7
sigma_V_p = 0 #0.0437
sigma_V_ap = 0 #0.0257
mean_tr = 0 #6.5e-11
sigma_tr = 0 #15.5E-12
mean_tf = 0 #6.57E-11
sigma_tf = 0 #4.66E-12
#the voltage sweep function is unabled 
if(initial_state):
	voltage_start = V_p
else:
	voltage_start = V_ap
voltage_end = -0.75
voltage_step = -0.1



fs = open (finished_results,'r')
content = fs.readlines()
fs.close()
exist_voltage= list()
exist_pulse = list()
exist_probability = list()
exist_trials = list()
for line in content:
	nums = line.split()
	exist_voltage.append( nums[0])
	exist_pulse.append(nums[1])
	exist_probability.append(nums[2])
	exist_results.append( "%.2f" % (float(nums[0])*10000000 + float(nums[1])) )
	exist_trials.append( float(nums[3]) )
for i_voltage in range(int( (voltage_end - voltage_start)/voltage_step) +1 ):
    for i_pulse in range( int( (pulse_end - pulse_start)/pulse_step) +1):
	voltage = voltage_start + voltage_step * i_voltage
	pulse = pulse_start + pulse_step * i_pulse

	fw1 = open(pulse_filename,'w')
	fw1.write("%.5g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g" % (voltage,pulse*1e-9,voltage,voltage,sigma_V_p,sigma_V_ap,mean_tr,sigma_tr,mean_tf,sigma_tf,sense_time,temperature));
	fw1.close()
	if( pulse >= pulse_start and pulse <= pulse_end ):
	  if ( ("%.2f" %  (voltage * 10000000 + pulse)) in exist_results ):
		ind = exist_results.index("%.2f" % (voltage * 10000000 + pulse))
		if ( float(exist_probability[ind]) > 1 - 100/float(trials) and float(trials) > exist_trials[ind] ):
			print ( str(voltage) + ' ' + str(pulse) +' is updating')
			cmd = './WERSim '+ str(trials) + ' ' + str(n_blocks) + ' ' + str(n_treads_p_block) + ' ' + str(initial_state) + ' '+pulse_filename + ' '+ str(STTorPS) + ' ./VInew/v_i_RTD_0.47u.txt  ./VInew/v_i_mos.txt 0 25e-15'
			print(cmd)
			os.system(cmd)
			if(initial_state == 0):
				fs = open ('p2ap.txt','r')
			else:
				fs = open ('ap2p.txt','r')
			allresult = fs.readlines()
			fs.close()
			thisresult = allresult[ len(allresult) -1]
			exist_probability[ind] = thisresult.split()[2]
			exist_trials[ind] = trials
			fw = open ( finished_results,'w')
			for i in range (len(exist_voltage)):
				fw.write( exist_voltage[i] + ' ' + exist_pulse[i] + ' '+ exist_probability[i]+' '+str(exist_trials[i])+'\n' )
			fw.close()
		else:	
			print ( str(voltage) + ' ' + str(pulse) +' '+ exist_probability[ind] +' in the '+ finished_results)
	  else:
		cmd = './WERSim '+ str(trials) + ' ' + str(n_blocks) + ' ' + str(n_treads_p_block) + ' ' + str(initial_state) + ' '+pulse_filename + ' ' + str(STTorPS) + '  ./VInew/v_i_RTD_0.47u.txt  ./VInew/v_i_mos.txt 0 25e-15'
		print(cmd)
		os.system(cmd)
		if(initial_state == 0):
			fs = open ('p2ap.txt','r')
		else:
			fs = open ('ap2p.txt','r')
		allresult = fs.readlines()
		fs.close()
		thisresult = allresult[ len(allresult) -1]
		probability = thisresult.split()[2]
		fw = open ( finished_results,'a')
		fw.write( str(voltage) + ' ' + str(pulse) + ' '+ probability+' '+str(trials)+'\n' )
		fw.close()
