#command [percent of slacks] [heterdesign, e.g. _GeSi without _
import os
import math
import string
import struct
import sys 

#usage python this [ndr_File] [read_config] [trials]
args = sys.argv
if(len(args)<4):
	print "usage: $python $0 [VI of ndr] [read_config] [# of trials]"
	quit()
ndr_file = args[1]
tempchar = ndr_file.replace('u','').split('.txt')[0].split('_')
ndr_width = float(tempchar[len(tempchar)-1])
fr = open(args[2],'r')
read_config = fr.readlines()
fr.close()
base_dir = './'
os.chdir(base_dir)
trials = int(args[3])
n_blocks = 128
n_treads_p_block = 128
initial_state = 1 # 1: ap2p, 0: p2ap
ndr_mode = 3 #0: 2: normal read, 3: ndr read

STTorPS = 0 # 0: STT, 1:Precessional Swiching

pulse_start = 0.1
pulse_end = 0.2
pulse_step = 0.1
exist_results = list()
finished_results = './finished_results.txt'
if (not os.path.isfile(finished_results)):
	fw = open(finished_results,'w')
	fw.write("ndr_width pulse_rise_time probability #_trials sensing_time Bitline_load Bitline_voltage std_Bitline_voltage\n")  
	fw.close()
template_pulse_filename = './template_pulse.txt'
pulse_filename = './var_input.txt'
temperature = 300;
V_p = -0.6
V_ap = -0.6
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
voltage_end = -0.7
voltage_step = -0.1



fs = open (finished_results,'r')
content = fs.readlines()
fs.close()
if(len(content)<1):
	fw.write("ndr_width pulse_rise_time probability #_trials sensing_time Bitline_load Bitline_voltage std_Bitline_voltage\n")
	
exist_voltage= list()
exist_pulse = list()
exist_probability = list()
exist_trials = list()
exist_sense_time = list()
exist_cload = list()
exist_Vbline = list()
exist_stdVndr = list()
for line in content[1:]:
	nums = line.split()
	exist_voltage.append( nums[0])
	exist_pulse.append(nums[1])
	exist_probability.append(nums[2])
	exist_results.append( "%.2f" % (float(nums[0])*1e7 + float(nums[1])*1e5 + float(nums[4])*1e11 + float(nums[5])*1e18 ) )
	exist_trials.append( float(nums[3]) )
	exist_sense_time.append( nums[4])
	exist_cload.append( nums[5])
	exist_Vbline.append(nums[6])
	exist_stdVndr.append(nums[7])
for line in read_config:
    parts = line.rstrip().split()
    voltage = float(parts[0])
    cload = float(parts[2])
    sense_time = float(parts[1])
    if(ndr_mode == 3): #ndr read
		mean_tr = cload/20e-15 * 0.1e-9
    if(ndr_mode == 2): #normal read
		mean_tr = cload/20e-15 * 0.1e-9
    for i_pulse in range( int( (pulse_end - pulse_start)/pulse_step) +1):
		pulse = pulse_start + pulse_step * i_pulse + mean_tr*1e9
		fw1 = open(pulse_filename,'w')
		fw1.write("%.5g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g\n%.4g" % (voltage,pulse*1e-9,voltage,voltage,sigma_V_p,sigma_V_ap,mean_tr,sigma_tr,mean_tf,sigma_tf,sense_time,temperature));
		fw1.close()
    if ( ("%.2f" %  (ndr_width * 1e7 + pulse*1e5 + cload*1e18 + sense_time*1e11)) in exist_results ):
		ind = exist_results.index("%.2f" % (ndr_width * 1e7 + pulse*1e5 + cload*1e18 + sense_time*1e11))
		if ( ( float(exist_probability[ind]) > 1 - 100/float(trials) or float(exist_probability[ind]) < 100/float(trials)) and float(trials) > exist_trials[ind] ):
			print ( str(ndr_width) + ' ' + str(pulse) +' is updating')
			cmd = './WERSim '+ str(trials) + ' ' + str(n_blocks) + ' ' + str(n_treads_p_block) + ' ' + str(initial_state) + ' '+pulse_filename + ' '+ str(STTorPS) + ' ' + ndr_file + ' v_i_mos.txt '+ str(ndr_mode) + ' ' + str(cload) + ' 2>&1 | tee sim.log'
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
			fs = open('sim.log','r')
			tempoutput = fs.readlines()
			fs.close()
			ave_endVndr = '-1'
			std_Vndr = '-1'
			for line in tempoutput:
				if('Average sensing margin is' in line):
					parts = line.rstrip().split()
					ave_endVndr = parts[len(parts)-1]
				if('Standard deviation of sensing margin is:' in line):
					parts = line.rstrip().split()
					std_Vndr = parts[len(parts)-1]
			exist_Vbline[ind] = ave_endVndr
			exist_stdVndr[ind] = std_Vndr
			fw = open ( finished_results,'w')
			fw.write("ndr_width pulse_rise_time probability #_trials sensing_time Bitline_load Bitline_voltage std_Bitline_voltage\n")
			for i in range (len(exist_voltage)):
				fw.write( exist_voltage[i] + ' ' + exist_pulse[i] + ' '+ exist_probability[i]+' '+str(exist_trials[i])+' ' + exist_sense_time[i]+' '+exist_cload[i]+' '+exist_Vbline[i]+' '+exist_stdVndr[i]+'\n' )
			fw.close()
		else:	
			print ( str(ndr_width) + ' ' + str(pulse) +' '+ exist_probability[ind] +' in the '+ finished_results)
    else:
		cmd = './WERSim '+ str(trials) + ' ' + str(n_blocks) + ' ' + str(n_treads_p_block) + ' ' + str(initial_state) + ' '+pulse_filename + ' ' + str(STTorPS) + ' ' + ndr_file+' v_i_mos.txt '+str(ndr_mode)+ ' ' + str(cload) + ' 2>&1 | tee sim.log'
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
		fs = open('sim.log','r')
		tempoutput = fs.readlines()
		fs.close()
		ave_endVndr = '-1'
		std_Vndr = '-1'
		for line in tempoutput:
			if('Average sensing margin is' in line):
				parts = line.rstrip().split()
				ave_endVndr = parts[len(parts)-1]
			if('Standard deviation of sensing margin is:' in line):
				parts = line.rstrip().split()
				std_Vndr = parts[len(parts)-1]
			
		fw = open ( finished_results,'a')
		fw.write( str(ndr_width) + ' ' + str(pulse) + ' '+ probability+' '+str(trials)+ ' ' + str(sense_time) + ' ' + str(cload) + ' ' +ave_endVndr+' '+std_Vndr+'\n' )
		fw.close()
		exist_voltage.append( str(ndr_width))
		exist_pulse.append(str(pulse))
		exist_probability.append(probability)
		exist_results.append( "%.2f" % (ndr_width*1e7 + pulse*1e5 + sense_time*1e11 + cload*1e18 ) )
		exist_trials.append( trials )
		exist_sense_time.append( str(sense_time))
		exist_cload.append( str(cload))
		exist_Vbline.append(ave_endVndr)
