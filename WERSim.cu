// -------------------------------------------
// Precessional Switching Old Parameters
// Updated November 21 - 2013
// Vector coordinates: [x; y; z]
// -------------------------------------------
#include <iostream>
#include <math.h>
#include "LLG.cu"
#include <vector>
#include <fstream>
#include <iomanip>
using namespace std;
// -------------------------------------------
// Calculate Wall time
// -------------------------------------------


//#define DEBUG

#include <sys/time.h>
#include "./Demagnetization_factors.cu"
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}



// -------------------------------------------
// WER Calculation
// -------------------------------------------
static cudaError_t crc;

void g_allocate_1D(double **g_f, int nsize, int *irc) {
/* allocate global double memory on GPU, return pointer to C */
   void *gptr;
   crc = cudaMalloc(&gptr,sizeof(double)*nsize);
   if (crc) {
      printf("cudaMalloc double Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *g_f = (double *)gptr;
   return;
}

void g_allocate_2D(double ***g_f, size_t * pitch, int row, int column, int *irc) {
/* allocate global double memory on GPU, return pointer to C */
   void *gptr;
   crc = cudaMallocPitch( &gptr, pitch, sizeof(double)*column, row);
   if (crc) {
      printf("cudaMalloc double Error=%d:%s,row=%d,column=%d\n",crc,
              cudaGetErrorString(crc),row,column);
      *irc = 1;
   }
   *g_f = (double **)gptr;
   return;
}


void copyin_gmemptr_1D (double *f, double *g_f, int nsize) {
/* copy double array from main memory to global GPU memory */
   crc = cudaMemcpy((void *)g_f,f,sizeof(double)*nsize,
                    cudaMemcpyHostToDevice);
   if (crc) {
      printf("cudaMemcpyHostToDevice1D double Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

void copyin_gmemptr_2D (double **f, double **g_f, size_t& pitch, int row, int column) {
/* copy double array from main memory to global GPU memory */
   crc = cudaMemcpy2D(g_f,pitch,f,sizeof(double)*column,sizeof(double)*column, row, cudaMemcpyHostToDevice);
   if (crc) {
      printf("cudaMemcpyHostToDevice2D double Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

void copyout_gmemptr_1D(double *f, double *g_f, int nsize) {
/* copy double array from global GPU memory to main memory */
   crc = cudaMemcpy(f,g_f,sizeof(double)*nsize,
                    cudaMemcpyDeviceToHost);
   if (crc) {
      printf("cudaMemcpyDeviceToHost1D double Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

void copyout_gmemptr_2D(double **f, double **g_f, size_t &pitch, int row, int column) {
/* copy double array from main memory to global GPU memory */
   crc = cudaMemcpy2D(f,sizeof(double)*column,g_f,pitch,sizeof(double)*column, row, cudaMemcpyDeviceToHost);
   if (crc) {
      printf("cudaMemcpyDeviceToHost2D double Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

int main(int argc, char* argv[])
{
	if(argc < 11) {
		cout<<" arguments: [trials] [number of blocks (32n)] [number of threads in block (32m)]  [initial state 0 = parallel, 1 = anti-parallel] [input pulse shape file] [enable VCMA?, 1:enable, 0: disable] [v_i characteristics of NDR] [v_i characteristics of MOS] [ndr write/read? 0: ndr write without ndr, 1: ndr write, 2: normal read without ndr, 3: ndr read; rise time is required except option 1] [Cbitline]"<<endl;
	return 1;
	}
	//dimention 
        double length = 50e-9;                            // MTJ length [m]
	double width = 50e-9;                            // MTJ width [m]
	double tfl = 1.18e-9;//origin:1.1e-9                               // Free layer thickness [m]
	double sigma_l = 1e-9;
	double sigma_w = 1e-9;
	double sigma_tfl = 0.003e-9;
	double sigma_mgo = 0.003e-9;
	double* lin_dep_factor = new double[9];
	double Cload = atof(argv[10]);
	//Demagnetization calculation
	cout<<"Start Demagnitization calculation"<<endl;
	double Nx = 0, Ny = 0, Nz = 0;
	Extract_linear_dependent(length, sigma_l,  width, sigma_w, tfl, sigma_tfl, Nx, Ny, Nz,lin_dep_factor);

#ifdef DEBUG
	std::cout<<Nx<<" "<<Ny<<" "<<Nz<<" "<<endl;
#endif
		
	double * g_lin_dep_factor;

	int irc = 0;
	g_allocate_1D(&g_lin_dep_factor, 9 ,&irc);
	if(irc!=0){
		cout<<"error in allocating memory in GPU"<<endl;
		return 1;
	}
	copyin_gmemptr_1D(  lin_dep_factor , g_lin_dep_factor, 9);
	
	fstream fs,fi;
	fi.open(argv[5],std::fstream::in);
	vector<double> v_par;
	double readBuffer,voltage=0;
	while(fi >> readBuffer){
		v_par.push_back(readBuffer);
	}
	fi.close();
	voltage = v_par[0]; 
	//Move pulse shape to 1D array to copy to GPU
	double *v_para;
	unsigned int n_par = v_par.size();
	v_para = new double[n_par];
	for (int i_p = 0; i_p < n_par; i_p ++ ){
		v_para[i_p] = v_par[i_p];
	}
	double t_pulse = v_par[1];
	vector<double>().swap(v_par);
	
	double t_step = 3e-12;                             // Time step [s]
	
	int trials = atoi(argv[1]);
	int GridSize = atoi(argv[2]);
	int GridLength = sqrt(GridSize);
	int BlockSize = atoi(argv[3]);
	bool enableVCMA=true; // is it precessional switching or STT
	if ( atoi(argv[6]) == 0){
		 enableVCMA = false;
	}
	
	int BlockLength = sqrt(BlockSize);
	int trials_p_thread = ceil(double(trials)/ double(GridSize*BlockSize));
	int real_trials = trials_p_thread * GridSize*BlockSize;
	cout<< "real number of trials is: "<<real_trials<<endl;
	dim3 dimBlock(BlockSize, 1);
	dim3 dimGrid(GridSize,1);
	int initial_state = atoi(argv[4]);
	if( initial_state == 1){
		fs.open("ap2p.txt",std::fstream::out | std::fstream::app);
	}
	else{
		fs.open("p2ap.txt", std::fstream::out | std::fstream::app);
	}
//Setup random variable

	
	double start_t = get_wall_time();
	cout<<"Start copy from host memory to GPU..."<<endl;

	double* global_writeSuccess = new double[real_trials];
	double * g_writeSuccess, *g_v_para;
	g_allocate_1D(&g_writeSuccess, real_trials ,&irc);
	copyin_gmemptr_1D( global_writeSuccess ,g_writeSuccess, real_trials);
	g_allocate_1D(&g_v_para, n_par, &irc);
	if(irc!=0){
		cout<<"error in allocating memory in GPU"<<endl;
		return 1;
	}
	copyin_gmemptr_1D( v_para, g_v_para, n_par);

/***************************************/
//Edition for NDR starts here	
	double* global_Energy = new double[real_trials];
	double* global_SwitchingTime = new double[real_trials];
	double* global_EndVndr = new double[real_trials];
	double *g_Energy, *g_SwitchingTime, *g_EndVndr, *VIGndr, *VIGmos, *g_VIGndr, *g_VIGmos; // the v and I characteristics of NDR and MOSFET
	int isNDR = atoi(argv[9]); // is NDR is added to the circuit
	int npoint_ndr = N_point; 
	VIGndr = Read_voltage_current(argv[7]); // Read in v_i characteristics
	VIGmos = Read_voltage_current(argv[8]);
//	size_t this_pitch;
//	g_allocate_2D(&g_VIGndr, &this_pitch, npoint_ndr, 3, &irc);
	g_allocate_1D(&g_VIGndr, npoint_ndr*3, &irc);
	g_allocate_1D(&g_VIGmos, npoint_ndr*3, &irc);
	if(irc!=0){
		cout<<"error in allocating memory in GPU"<<endl;
		return 1;
	}
	copyin_gmemptr_1D( VIGndr, g_VIGndr, npoint_ndr*3);
	copyin_gmemptr_1D( VIGmos, g_VIGmos, npoint_ndr*3);
	g_allocate_1D(&g_Energy, real_trials ,&irc);
	g_allocate_1D(&g_SwitchingTime, real_trials ,&irc);
	g_allocate_1D(&g_EndVndr, real_trials ,&irc);
	if(irc!=0){
		cout<<"error in allocating memory in GPU"<<endl;
		return 1;
	}
	copyin_gmemptr_1D( global_Energy ,g_Energy, real_trials);
	copyin_gmemptr_1D( global_SwitchingTime ,g_SwitchingTime, real_trials);
	copyin_gmemptr_1D( global_EndVndr ,g_EndVndr, real_trials);

#ifdef OUTPUT_DETAIL	
	double* global_NDRturnR = new double[real_trials];
	double* global_initialR = new double[real_trials];
	double* global_NDRoffR = new double[real_trials];
	double *g_NDRturnR, *g_initialR, *g_NDRoffR ;
	g_allocate_1D(&g_NDRturnR, real_trials ,&irc);
	g_allocate_1D(&g_initialR, real_trials ,&irc);
	g_allocate_1D(&g_NDRoffR, real_trials ,&irc);
	if(irc!=0){
		cout<<"error in allocating memory in GPU"<<endl;
		return 1;
	}
	copyin_gmemptr_1D( global_NDRturnR ,g_NDRturnR, real_trials);
	copyin_gmemptr_1D( global_initialR ,g_initialR, real_trials);
	copyin_gmemptr_1D( global_NDRoffR ,g_NDRoffR, real_trials);
#endif
// Editional for NDR ends here	
/***************************************/

	double copyin_t = get_wall_time();
	cout<<"Start GPU calculation..."<<endl;
#ifdef OUTPUT_DETAIL
   	LLG<<<dimGrid,dimBlock>>>(g_v_para, g_writeSuccess, initial_state, t_step, trials_p_thread, enableVCMA, length, sigma_l,  width, sigma_w, tfl, sigma_tfl, sigma_mgo, Nx, Ny, Nz, g_lin_dep_factor,isNDR,g_VIGndr,g_VIGmos,  g_Energy, g_SwitchingTime, g_EndVndr, Cload, g_NDRturnR, g_initialR,g_NDRoffR);
#endif
#ifndef OUTPUT_DETAIL
   	LLG<<<dimGrid,dimBlock>>>(g_v_para, g_writeSuccess, initial_state, t_step, trials_p_thread, enableVCMA, length, sigma_l,  width, sigma_w, tfl, sigma_tfl, sigma_mgo, Nx, Ny, Nz, g_lin_dep_factor,isNDR,g_VIGndr,g_VIGmos,  g_Energy, g_SwitchingTime, g_EndVndr,Cload);
#endif
	cudaDeviceSynchronize();
	double calculate_t = get_wall_time();
	cout<<"Start copy GPU to host memory..."<<endl;
	copyout_gmemptr_1D(global_writeSuccess, g_writeSuccess, real_trials);
/****************************/
//Edition for NDR starts here
	copyout_gmemptr_1D(global_Energy, g_Energy, real_trials);
	copyout_gmemptr_1D(global_SwitchingTime, g_SwitchingTime, real_trials);
	copyout_gmemptr_1D(global_EndVndr, g_EndVndr, real_trials);
	cout<<"Start couting final result...\n**********************"<<endl;
#ifdef OUTPUT_DETAIL
	copyout_gmemptr_1D(global_NDRturnR, g_NDRturnR, real_trials);
	copyout_gmemptr_1D(global_initialR, g_initialR, real_trials);
	copyout_gmemptr_1D(global_NDRoffR, g_NDRoffR, real_trials);
	fstream fout;
	fout.open("output_detailR.txt",std::fstream::out);
	for (int i_trial = 0; i_trial<real_trials ; i_trial++){
		fout << "switched/R_initial/NDRturn/NDRoff\t"<<global_writeSuccess[i_trial]<<"\t"<<global_initialR[i_trial]<<"\t"<<global_NDRturnR[i_trial]<<"\t"<<global_NDRoffR[i_trial]<<std::endl;
	}
	fout.close();
	
#endif
	double total_energy = 0;
	double total_switchingtime = 0;
	double ave_vndr = 0;
	int read_failures = 0;
	double min_margin = 0;
	if(n_par > 12) {
		min_margin = v_para[12];
	}
	for(int i_e = 0; i_e < real_trials; i_e++){
		total_energy += global_Energy[i_e];
		total_switchingtime += global_SwitchingTime[i_e];
		ave_vndr += global_EndVndr[i_e]; // Sensed voltage
		read_failures += (global_EndVndr[i_e] > min_margin)? 0:1;
	}
	ave_vndr /= real_trials;
	if(isNDR >=2){ // in read mode
		cout<<"Results of read mode ("<<isNDR<<")"<<endl;
		double std_vndr = 0;
		for( int i_trial = 0; i_trial < real_trials; i_trial++){
			std_vndr += (global_EndVndr[i_trial] - ave_vndr)*(global_EndVndr[i_trial] - ave_vndr);
		}
		std_vndr = sqrt(std_vndr/real_trials);
		cout<<"  Average sensing margin is: "<<ave_vndr<<endl;
		cout<<"  Standard deviation of sensing margin is: " << std_vndr<<endl;
		cout<<"  Number of read failures (< minimum margin of "<<min_margin<<"): "<<read_failures<<endl;
	}else{
		cout<<"Results of write mode ("<<isNDR<<")"<<endl;
		cout<<"  Average switching power is: "<< total_energy/real_trials <<"\nAverage switching time is: "<<total_switchingtime/real_trials<<endl;
		cout<<"  Average Vndr after switching is: "<<ave_vndr<<endl;
	}

//Edition for NDR ends here
/***************************/
	double copyout_t = get_wall_time();

	int sum=0;
   	for( int i =0 ; i< real_trials; i++){
   		sum+= global_writeSuccess[i];
   	}
	double couting_t = get_wall_time();
	cout << "**********************\nMonte-Carlo Simulation Results:\n  switching from "<< ((initial_state==0)? "P" : "AP")<<" to "<< ((initial_state==1)? "P" : "AP")<<"\n"
	<<"  Pulse voltage: "<<voltage<<" V, time step: "<<t_step<<"s\n  "
        <<sum<<" trials success out of total "<<real_trials<<" trials, switching rate: "<<std::setprecision(9)<<(double(sum) / double(real_trials)) << endl;
	cout<<"**********************\nRuntime summary:\n"<<" copy memory to GPU: "<<copyin_t - start_t
<<" s, GPU calculation: "<<calculate_t - copyin_t<<" s, copy memory out to CPU: "
<<copyout_t - calculate_t <<" s, couting switching: "<< couting_t - copyout_t <<" s."<<endl;
	fs << voltage <<" "<<t_pulse<<" "<<std::setprecision(9)<<(double(sum) / double(real_trials))<<" "<<couting_t-start_t<<" s"<<endl;
	fs.close();
	return 0;
}





