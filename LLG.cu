#include <stdio.h>      
#include <math.h>    
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include </u/local/cuda/5.0/include/cuda.h>
#include </u/local/cuda/5.0/include/cuda_runtime.h>
#include </u/local/cuda/5.0/include/curand_kernel.h>
#include "./NDR_Solver.cu"
using namespace std;
#define VARIATION
//#define OUTPUT_DETAIL
#define shaodi_pi 3.1415926
#define CUDA_CALL(x) do { if( (x) ! =  cudaSuccess ){\
	printf("Error at %s:%d\n",__FILE__,__LINE__ );\
	exit(1);} } while(0) 
#define CURAND_CALL(x) do { if((x)!=CURAND_STATUS_SUCCESS) {\
	printf("Error at %s:%d\n",__FILE__,__LINE__);\
	exit(1);}} while(0)

#ifndef OUTPUT_DETAIL
__global__ void LLG(double* g_v_para, double* writeSuccess, int initialstate, double t_step, 
int trials_p_thread, bool isPS, double ori_length, double sigma_l, double ori_width, double sigma_w, 
double ori_tfl, double sigma_tfl, double sigma_mgo, double ori_Nx, double ori_Ny, double ori_Nz, 
double* g_lin_dep_factor, int isNDR, const double* VIGndr, const double* VIGmos, double* g_Energy, 
double* g_SwitchingTime, double* g_EndVndr, double Cload ){
#endif
#ifdef OUTPUT_DETAIL
__global__ void LLG(double* g_v_para, double* writeSuccess, int initialstate, double t_step,
int trials_p_thread, bool isPS, double ori_length, double sigma_l, double ori_width, 
double sigma_w, double ori_tfl, double sigma_tfl, double sigma_mgo, double ori_Nx, 
double ori_Ny, double ori_Nz, double* g_lin_dep_factor, int isNDR, const double* VIGndr, 
const double* VIGmos, double* g_Energy, double* g_SwitchingTime, double* g_EndVndr, double Cload,
double* g_NDRturn, double* g_initialR, double* g_NDRoff){
#endif
/* -------------------------------------------
 Input Parameters From User
 -------------------------------------------*/
int this_id = (blockIdx.x * blockDim.x + threadIdx.x) ;
//initiate state for following random generation
curandState_t localState;
curand_init(this_id, this_id, 0, &localState);

double Nx = ori_Nx;//origin:                                                // x Demagnetization factor
double Ny = ori_Ny;                                                // y Demagnetization factor
double Nz = ori_Nz;                                                // z Demagnetization factor

//Parameter calculation
double length = ori_length; //length of MTJ
double width = ori_width; //width of MTJ
double Rp0 = 2e3;
double dMgO_a = 1.54e-3, dMgO_b = 1.1537278e10;//origin:9.24e9;
double Area = shaodi_pi*length*width/4;            // Area without variation
double dMgO_base = (log(Rp0 * Area * 10e12) - log(dMgO_a)) / dMgO_b;     // MgO thickness [m]
double TMR = 1.5;                                    // TMR at zero bias voltage
double Rap0 = Rp0 *( 1+TMR);
double Temperature = g_v_para[11];//27+273;                       // Temperature
double pulse_width = g_v_para[1];
double V_p = g_v_para[2], V_ap = g_v_para[3];
double sigma_V_p = g_v_para[4], sigma_V_ap = g_v_para[5];
if( sigma_V_p == 0){
	sigma_V_p = 1e-9;
}
double mean_tr = g_v_para[6], sigma_tr = g_v_para[7], mean_tf = g_v_para[8], sigma_tf= g_v_para[9], delay_time = 0e-9, sense_time = g_v_para[10]; 

int n_sim = (pulse_width+delay_time + sense_time)/t_step ;//Simulation time

//double initial_Vndr = 0;
#ifdef OUTPUT_DETAIL
double	peak_voltage = Peak_voltage(VIGndr); 
double	peak_current = IG_V(peak_voltage,VIGndr,1);
#endif

/* -------------------------------------------
 Constants
 -------------------------------------------*/

double hbar = 1.05457173e-34;                                      // Reduced Planck constant, [J*s]
double k = 1.3806488e-23;                                          // Boltzmann constant, [J/K]
double u0 = 4e-7*shaodi_pi;                                               // Vacuum permeability, [V�s/(A�m)]
double q = 1.60217657e-19;                                         // Electron charge, [C]
double alphac = 0.02;                                              // LLGE damping factor
double gammap = (221276/(1+pow(alphac,2)));                             // Gyromagnetic ratio [m/(A x s)]
double T0 = 1120;
double Ms0 = 1393128.323;//origin:1.44e6;
double Ki0 =1.479036e-3;//origin:1.46e-3;
double Xi0 = 0; //53.39247e-15; //origin:58.9e-15;
if(isPS) Xi0 = 53.39247e-15;
double P_tunnel = 0.2;                                  // the polarization of the tunnel currentdouble
double Pol = 0.6;                                                  // Polarization for Spin Torque

/******************simulation trials *************/
for( int i_trial = 0; i_trial < trials_p_thread; i_trial++){

length = ori_length; //length of MTJ
width = ori_width; //width of MTJ
double tfl = ori_tfl; //thickness of free layer

double rise_time = mean_tr;
double fall_time = mean_tf;
#ifdef VARIATION
//Dimention variation
rise_time += sigma_tr*curand_normal_double(&localState);
fall_time += sigma_tf*curand_normal_double(&localState);
double v_variation = sigma_V_p * curand_normal_double(&localState);
length = ori_length + sigma_l*curand_normal_double(&localState); 
width = ori_width + sigma_w*curand_normal_double(&localState);
tfl = ori_tfl + sigma_tfl*curand_normal_double(&localState);
double dMgO = dMgO_base + sigma_mgo*curand_normal_double(&localState);
double temp_Nx = ori_Nx + g_lin_dep_factor[0]*(length-ori_length) + g_lin_dep_factor[1] * ( width - ori_width) + g_lin_dep_factor[2] * (tfl - ori_tfl) ;
double temp_Ny = ori_Ny + g_lin_dep_factor[3]*(length-ori_length) + g_lin_dep_factor[4] * ( width - ori_width) + g_lin_dep_factor[5] * (tfl - ori_tfl) ;
double temp_Nz = ori_Nz + g_lin_dep_factor[6]*(length-ori_length) + g_lin_dep_factor[7] * ( width - ori_width) + g_lin_dep_factor[8] * (tfl - ori_tfl) ;
Nx = temp_Nx / (temp_Nx + temp_Ny + temp_Nz);
Ny = temp_Ny / (temp_Nx + temp_Ny + temp_Nz);
Nz = temp_Nz / (temp_Nx + temp_Ny + temp_Nz); 
#endif

Area = shaodi_pi*length*width/4;            // Area without variation
double areamtj = Area  ;                                // MTJ area [m^2]
double Rp = exp(dMgO * dMgO_b)*dMgO_a / (Area * 10e12);
double Rap = (1+TMR)*Rp;                                           // Anti-parallel resistance [Ohms]
double B1 = 0;//origin 0.2                                   // Field-like torque linear parameter [unitless]
double B2 = 0;//origin 0.02;                                     // Field-like torque quadratic parameter [1/A]
int initial_state = initialstate;                          // Inital state [0 = parallel, 1 = anti-parallel]
double P [3] = {0, 0, -1};                             // Direction of polarization
double Ext [3] = {0, 0, 0};                        // External magnetic field [A/m] - 1 oersted [Oe] = 79.5774715459424 ampere/meter [A/m]


//double t_delay = 2e-9;                             // Time to initiate pulse application [s]

double Ms = Ms0 * ( 1 - pow(Temperature/T0,1.5));                  // Saturation magnetization [A/m] - 1e6 A/m = 1000 emu/cc 
double dstray = 20e-9, tstray = 1.164656e-9;
//double Ext[3]	 = {-Ms*length*width/4/shaodi_pi*((dstray+tstray)/(pow(length/2,2)*sqrt(pow(length/2,2)+pow(dstray+tstray,2)))-(dstray-tstray)/(pow(length/2,2)*sqrt(pow(length/2,2)+pow(dstray-tstray,2)))),0,0};
double Ki = Ki0 * pow(Ms/Ms0, 2.18);                          // Anisotropy field constant [J/m^2]
double Xi = Xi0* pow(Ms/Ms0, 2.83);                                // VCMA field constant [J/(V x m)]
double Gt = 1/(Rp*(1+(TMR/(TMR+2))));                              // Direct elastic tunneling conductance [S]
double KiPF = (2*Ki)/(tfl*u0*Ms);                                  // Prefactor for interface anisotropy effective field
double VCMAPF = (2*Xi)/(u0*Ms*dMgO*tfl);                           // Prefactor for VCMA effective field
double Gsi	= 0;                                                    // Conductance due to imperfections in Mgo [S]

//double Jc0 = (2*Ms*tfl*q*u0)/(hbar*Pol);                           // Normalization Constant for Current Density

double volume = areamtj*tfl;                                       // MTJ volume [m^3]
double Hth = sqrt((2*k*Temperature*alphac)/(u0*gammap*Ms*volume*t_step));    // Amplitude of Thermal Field
//int this_id = (blockIdx.x * blockDim.x + threadIdx.x) * trials_p_thread + i_trial;

/* -------------------------------------------
 Internal Variables
 -------------------------------------------*/

double costheta = 0;                                       // the angle between the magnization  of free and reference layers
double g_sv = 0;                                        // the polarization efficiency in spin valve
double g_tunnel = 0;                                    // the polarization efficiency in tunnel current


//double m_old [3] = {0, 0, 0};                              // Normalized previous magnetization
double Heff_old [3] = {0, 0, 0};                           // Previous Heff components [A/m]
double m_int [3] = {0, 0, 0};                              // Intermediate normalized magnetization
double dm_int [3] = {0, 0, 0};                             // Intermediate derivative of magnetization
double M_int [3] = {0, 0, 0};                              // Intermediate denormalized magnetization                              
//double Heff_int [3] = {0, 0, 0};                           // Intermediate Heff components [A/m]
double dm [3] = {0, 0, 0};                                 // Time derivative of magnetization [1/s]
double M [3] = {0, 0, 0};                                  // Denormalized magnetization
double mcrossp_int [3] = {0, 0, 0};                        // Intermediate cross product components (m x p)
double mcrossHeff_int [3] = {0, 0, 0};                     // Intermediate cross product components (m x Heff)
double mcrossHth_int [3] = {0, 0, 0};                      // Intermediate cross product components (m x Hth)
double mcrossmcrossp_int [3] = {0, 0, 0};                  // Intermediate double cross product components (m x m x p)
double mcrossmcrossHeff_int [3] = {0, 0, 0};               // Intermediate double cross product components (m x m x Heff)
double mcrossp [3] = {0, 0, 0};                            // Cross product components (m x p)
double mcrossHeff [3] = {0, 0, 0};                         // Cross product components (m x Heff)
double mcrossHth [3] = {0, 0, 0};                          // Cross product components (m x Hth)
double mcrossmcrossp [3] = {0, 0, 0};                      // Cross product components (m x m x p)
double mcrossmcrossHeff [3] = {0, 0, 0};                   // Cross product components (m x m x Heff)
double randomHth [3] = {0, 0, 0};                          // Vector of random variables
double STT  = 0;                                        // Strenght of STT term
double FLT  = 0;                                        // Strenght of FLT term

// -------------------------------------------
// Initialize Variables
// -------------------------------------------
double m [3] = {0, 0, 1};                             // Normalized mangetization
double R  = Rap;                						// MTJ resistance [Ohms]
if(initial_state != 1){
                       
    R  = Rp;                                    // MTJ resistance [Ohms]
    m[2]  = -1;                             // Normalized mangetization
}

double J = 0;                                          // Current density [A/m^2]
double V = 0;                                          // MTJ Voltage [V]

 
double V_offset = 0;//1e-10 * (this_id * trials_p_thread+i_trial);

/*********** edition for NDR starts here ***********/
//The parameters for calculating NDR
double Vndr = 0; 
double Imtj = 0; //current through MTJ and nmos
double Vmos = 0; 
double d_Rmtj = 0; // delta Rmtj
double d_Imtj = 0; // delta Imtj
double Indr = 0; //current through NDR
double d_Vndr = 0; //
double Csline = Cload;
double Cbline = Cload;
double d_vdd =0, new_vdd =0, vdd =0;
V = 0; // mtj voltage
#ifdef OUTPUT_DETAIL
	bool isNDRturn = false, isNDRoff = false;
	g_NDRturn [this_id*trials_p_thread + i_trial] = 0;
	g_NDRoff [this_id*trials_p_thread + i_trial] = 0;
	g_initialR [this_id*trials_p_thread + i_trial] = R;
#endif
if(isNDR == 1){ //NDR write
	//Vndr = initial_Vndr;
	Vndr = Solve_stable_vndr(VIGndr, VIGmos, R, V_ap);
	Indr = IG_V(Vndr, VIGndr,1);
	Imtj = Indr;
	Vmos = V_I(Imtj,VIGmos);
	vdd = V_ap;
}
double energy = V_ap*V_ap*Cload; // Pre-charge energy
bool isSwitched = false;
g_SwitchingTime[this_id * trials_p_thread+i_trial] = pulse_width;
/*********** edition for NDR ends here ***********/

for(int i=1;i<=n_sim;i++){

    // Update values
    double m_old [3] = {m[0], m[1], m[2]};
     
    
    // Update voltage/current density
    double V_ub = V_ap + V_offset;
#ifdef VARIATION
    V_ub += v_variation * ( 1 + (sigma_V_ap/sigma_V_p - 1) * (R - Rp0)/(Rap0 - Rp0) ) ;
#endif
    double curr_time = i * t_step;
    if(curr_time < delay_time || curr_time > delay_time + pulse_width){
	new_vdd = 0;
    }
    else{
	if(curr_time < delay_time + rise_time){
		new_vdd = (curr_time - delay_time)/rise_time * V_ub;
	}
	else{
	    if(curr_time <= delay_time + pulse_width - fall_time){
		new_vdd = V_ub;
	    }
	    else{
		new_vdd = V_ub * (  delay_time + pulse_width - curr_time) / fall_time;
	    }
	}
    }
    d_vdd = new_vdd - vdd;
    vdd = new_vdd;
//NDR calculation
/*********** edition for NDR starts here ***********/

#ifdef OUTPUT_DETAIL
    if(isNDR==1){
	if( !isNDRturn && abs(Vndr) > abs(peak_voltage)){
		g_NDRturn [this_id*trials_p_thread + i_trial] = R;
		isNDRturn = true;
	}
	if( isNDRturn && !isNDRoff && abs(Indr) < 0.25*abs(peak_current)){
		g_NDRoff [this_id*trials_p_thread + i_trial] = R;
		isNDRoff = true;
	}
    } else if(isNDR>=2){
	if( !isNDRturn && abs(Vndr) > abs(peak_voltage)){
		g_NDRturn [this_id*trials_p_thread + i_trial] = curr_time;
		isNDRturn = true;
	}
	if (!isNDRoff && curr_time >= delay_time + pulse_width+sense_time){
//	if (!isNDRoff && i>=1){
		isNDRoff = true;
//		g_initialR [this_id*trials_p_thread + i_trial] = d_Vndr;
//		g_NDRturn [this_id*trials_p_thread + i_trial] = Imtj;
	    	g_NDRoff [this_id*trials_p_thread + i_trial] = Vndr + V + Vmos;
	}		
    }
#endif
    if(isNDR==1 ){ //NDR write for AP-MTJ
	//Solve the series of one NDR, one MTJ and one MOS with cap at MTJ
	d_Imtj = ( t_step*Indr - t_step*Imtj - Cload*Imtj*d_Rmtj) / ( Cload/IG_V(Vmos,VIGmos,2) + Cload*R);
	Imtj += d_Imtj;
	Vmos = V_I(Imtj, VIGmos);
	Vndr = vdd - Imtj*R - Vmos;
	Indr = IG_V(Vndr,VIGndr,1);
	V = Imtj*R;
    	energy += vdd * Indr * t_step;
    }
    else if(isNDR == 3){ // NDR read
	if(curr_time < delay_time+pulse_width){ // precharging
            //solve equation: Imtj = (vdd - Vmos(Imtj) - Vndr)/Rmtj = d(Csline*Vndr)/dt + Indr(Vndr)
            Vmos = V_I(Imtj,VIGmos);
            d_Imtj = (d_vdd - d_Vndr - d_Rmtj*Imtj)/(R+1/IG_V(Vmos,VIGmos,2));
            Imtj += d_Imtj;
            Indr = IG_V(Vndr,VIGndr,1);
            d_Vndr = ( Imtj - Indr ) * t_step / Csline;
            Vndr += d_Vndr;
            V = Imtj*R;
    	    energy += vdd * Imtj * t_step;
        }
        else{ // discharging state
            //Solve equation: - d ((Imtj*R + Vmos(Imtj) + Vndr)*Cbline )/dt = Imtj = d(Vndr*Csline)/dt + Indr(Vndr)

            d_Imtj = (-Imtj*t_step - Cbline*Imtj*d_Rmtj -d_Vndr*Cbline) / ( Cbline/IG_V(Vmos,VIGmos,2) + Cbline*R);
            Imtj += d_Imtj;
            Vmos = V_I(Imtj,VIGmos);
            d_Vndr = (Imtj - Indr) * t_step / Csline;
            Vndr += d_Vndr;
            Indr = IG_V(Vndr,VIGndr,1);
            V = Imtj*R;
        }
    }
    else if (isNDR == 2 || isNDR == 0){
	 if(curr_time < delay_time+pulse_width){ // precharging
            //solve equation: Imtj = (vdd - Vmos(Imtj) )/Rmtj = 0
            d_Imtj = (d_vdd-d_Rmtj*Imtj)/(R+1/IG_V(Vmos,VIGmos,2));
            Imtj += d_Imtj;
            V = Imtj*R;
            Vmos = vdd - V;
	    Indr = Imtj + R*Cbline*d_Imtj/t_step + Imtj*Cbline*d_Rmtj/t_step + Cbline/IG_V(Vmos,VIGmos,2)*d_Imtj/t_step;//The total I but not current of ndr, because there is no ndr
	    energy += vdd * Indr * t_step;
        }
        else{ // discharging state
            //Solve equation: - d ((Imtj*R + Vmos(Imtj) )*Cbline )/dt = Imtj 
    
            d_Imtj = (-Imtj*t_step - Cbline*Imtj*d_Rmtj ) / ( Cbline/IG_V(Vmos,VIGmos,2) + Cbline*R);
            Imtj += d_Imtj;
            Vmos = V_I(Imtj,VIGmos);
            V = Imtj*R;
        }
    }
//    else{
	//Solve the series of one MTJ and one MOS with cap at MTJ
//	if(curr_time < delay_time+pulse_width){
//	    d_Imtj = (d_vdd-d_Rmtj*Imtj)/(R+1/IG_V(Vmos,VIGmos,2));
//            Imtj += d_Imtj;
//            V = Imtj*R;
//            Vmos = vdd - V;
//	    Indr = Imtj + R*Cload*d_Imtj/t_step + Imtj*Cload*d_Rmtj/t_step + Cload/IG_V(Vmos,VIGmos,2)*d_Imtj/t_step;//The total I but not current of ndr, because there is no ndr
//	    energy += vdd * Indr * t_step;
//    }
//Test whether switched
    if(!isSwitched){
	if( (initial_state ==0 && R >= Rp*(1+TMR/2)) || ( initial_state ==1 && R <= Rp*(1+TMR/2)) ){
	    isSwitched = true;
	    g_SwitchingTime[this_id * trials_p_thread+i_trial] = curr_time;
	}
    }
            
/*********** edition for NDR ends here ***********/

    
	

    // Update effective magnetic field Heff_old
    Heff_old[0] = Ext[0]-Ms*Nx*m_old[0];
    Heff_old[1] = Ext[1]-Ms*Ny*m_old[1];
    Heff_old[2] = Ext[2]-Ms*Nz*m_old[2]+(KiPF*m_old[2]-VCMAPF*m_old[2]*V);
    
    //Calculate STT factor
    J = V/(R*areamtj);
//    costheta = m_old[0]*P[0] + m_old[1]*P[1] + m_old[2]*P[2];
//    g_tunnel = 1/2 * P_tunnel / ( 1 + pow(P_tunnel,2)*costheta);
//    g_sv = 1 / ( -4 + pow(( 1 / sqrt(Pol) + sqrt(Pol) ), 3) * (3 + costheta) / 4);
//    STT = gammap*J* hbar*(g_tunnel+g_sv)/(2*Ms*tfl*q*u0);
    STT = Pol*gammap*J* hbar/(2*Ms*tfl*q*u0);
    //STT = gammap*J/Jc0;
    FLT = STT*B1+STT*B2*areamtj*J;

    // Calculate m x Hth
    mcrossHth_int[0]=m_old[1]*randomHth[2]-m_old[2]*randomHth[1];
    mcrossHth_int[1]=m_old[2]*randomHth[0]-m_old[0]*randomHth[2];
    mcrossHth_int[2]=m_old[0]*randomHth[1]-m_old[1]*randomHth[0];

    // Calculate m x p and m x m x p
    mcrossp_int[0]=m_old[1]*P[2]-m_old[2]*P[1];
    mcrossp_int[1]=m_old[2]*P[0]-m_old[0]*P[2];
    mcrossp_int[2]=m_old[0]*P[1]-m_old[1]*P[0];
    mcrossmcrossp_int[0]=m_old[1]*mcrossp_int[2]-m_old[2]*mcrossp_int[1];
    mcrossmcrossp_int[1]=m_old[2]*mcrossp_int[0]-m_old[0]*mcrossp_int[2];
    mcrossmcrossp_int[2]=m_old[0]*mcrossp_int[1]-m_old[1]*mcrossp_int[0];

    // Calculate m x Heff and m x m x Heff
    mcrossHeff_int[0]=m_old[1]*Heff_old[2]-m_old[2]*Heff_old[1];
    mcrossHeff_int[1]=m_old[2]*Heff_old[0]-m_old[0]*Heff_old[2];
    mcrossHeff_int[2]=m_old[0]*Heff_old[1]-m_old[1]*Heff_old[0];
    mcrossmcrossHeff_int[0]=m_old[1]*mcrossHeff_int[2]-m_old[2]*mcrossHeff_int[1];
    mcrossmcrossHeff_int[1]=m_old[2]*mcrossHeff_int[0]-m_old[0]*mcrossHeff_int[2];
    mcrossmcrossHeff_int[2]=m_old[0]*mcrossHeff_int[1]-m_old[1]*mcrossHeff_int[0];
    // Use the LLG equation w/ Heun's Method to update the magnetization
    dm_int[0] = -gammap*(mcrossHeff_int[0]+mcrossHth_int[0]) - gammap*alphac*mcrossmcrossHeff_int[0] + STT*mcrossmcrossp_int[0] + FLT*mcrossp_int[0];
    dm_int[1] = -gammap*(mcrossHeff_int[1]+mcrossHth_int[1]) - gammap*alphac*mcrossmcrossHeff_int[1] + STT*mcrossmcrossp_int[1] + FLT*mcrossp_int[1];
    dm_int[2] = -gammap*(mcrossHeff_int[2]+mcrossHth_int[2]) - gammap*alphac*mcrossmcrossHeff_int[2] + STT*mcrossmcrossp_int[2] + FLT*mcrossp_int[2];
    M_int[0] = m_old[0] + (dm_int[0]*t_step);
    M_int[1] = m_old[1] + (dm_int[1]*t_step);
    M_int[2] = m_old[2] + (dm_int[2]*t_step);
    m_int[0] = M_int[0]/sqrt(M_int[0]*M_int[0]+M_int[1]*M_int[1]+M_int[2]*M_int[2]);
    m_int[1] = M_int[1]/sqrt(M_int[0]*M_int[0]+M_int[1]*M_int[1]+M_int[2]*M_int[2]);
    m_int[2] = M_int[2]/sqrt(M_int[0]*M_int[0]+M_int[1]*M_int[1]+M_int[2]*M_int[2]);

    // Update the thermal field and current values (time evolves)
    
	double2 gen_x12;
	double gen_x3;
	gen_x12 = curand_normal2_double(&localState);
      randomHth[0] = Hth*gen_x12.x;
      randomHth[1] = Hth*gen_x12.y;
	gen_x3 = curand_normal_double(&localState);
      randomHth[2] = Hth*gen_x3;
    
//STT calculation
//    costheta = m_int[0]*P[0] + m_int[1]*P[1] + m_int[2]*P[2];
//    g_tunnel = 1/2 * P_tunnel / ( 1 + pow(P_tunnel,2)*costheta);
//    g_sv = 1 / ( -4 + pow(( 1 / sqrt(Pol) + sqrt(Pol) ), 3) * (3 + costheta) / 4); 
//    STT = gammap*J* hbar*(g_tunnel+g_sv)/(2*Ms*tfl*q*u0);
    STT = Pol*gammap*J* hbar/(2*Ms*tfl*q*u0);

    //STT = gammap*J/Jc0;
    FLT = STT*B1+STT*B2*areamtj*J;

    // Update intermediate effective magnetic field Heff
    double Heff_int [3] = {Ext[0]-Ms*Nx*m_int[0], Ext[1]-Ms*Ny*m_int[1], Ext[2]-Ms*Nz*m_int[2]+(KiPF*m_int[2]-VCMAPF*m_int[2]*V)};

    // Calculate m x Hth
    mcrossHth[0]=m_int[1]*randomHth[2]-m_int[2]*randomHth[1];
    mcrossHth[1]=m_int[2]*randomHth[0]-m_int[0]*randomHth[2];
    mcrossHth[2]=m_int[0]*randomHth[1]-m_int[1]*randomHth[0];
    // Calculate m x p and m x m x p
    mcrossp[0]=m_int[1]*P[2]-m_int[2]*P[1];
    mcrossp[1]=m_int[2]*P[0]-m_int[0]*P[2];
    mcrossp[2]=m_int[0]*P[1]-m_int[1]*P[0];
    mcrossmcrossp[0]=m_int[1]*mcrossp[2]-m_int[2]*mcrossp[1];
    mcrossmcrossp[1]=m_int[2]*mcrossp[0]-m_int[0]*mcrossp[2];
    mcrossmcrossp[2]=m_int[0]*mcrossp[1]-m_int[1]*mcrossp[0];

    // Calculate m x Heff and m x m x Heff
    mcrossHeff[0]=m_int[1]*Heff_int[2]-m_int[2]*Heff_int[1];
    mcrossHeff[1]=m_int[2]*Heff_int[0]-m_int[0]*Heff_int[2];
    mcrossHeff[2]=m_int[0]*Heff_int[1]-m_int[1]*Heff_int[0];
    mcrossmcrossHeff[0]=m_int[1]*mcrossHeff[2]-m_int[2]*mcrossHeff[1];
    mcrossmcrossHeff[1]=m_int[2]*mcrossHeff[0]-m_int[0]*mcrossHeff[2];
    mcrossmcrossHeff[2]=m_int[0]*mcrossHeff[1]-m_int[1]*mcrossHeff[0];

    // Now use intermediate value in final value computation 
    dm[0] = -gammap*(mcrossHeff[0]+mcrossHth[0]) - gammap*alphac*mcrossmcrossHeff[0] + STT*mcrossmcrossp[0] + FLT*mcrossp[0];
    dm[1] = -gammap*(mcrossHeff[1]+mcrossHth[1]) - gammap*alphac*mcrossmcrossHeff[1] + STT*mcrossmcrossp[1] + FLT*mcrossp[1];
    dm[2] = -gammap*(mcrossHeff[2]+mcrossHth[2]) - gammap*alphac*mcrossmcrossHeff[2] + STT*mcrossmcrossp[2] + FLT*mcrossp[2];
    M[0] = m_old[0] + (t_step/2)*(dm[0] + dm_int[0]);
    M[1] = m_old[1] + (t_step/2)*(dm[1] + dm_int[1]);
    M[2] = m_old[2] + (t_step/2)*(dm[2] + dm_int[2]);
    m[0] = M[0]/sqrt(M[0]*M[0]+M[1]*M[1]+M[2]*M[2]);
    m[1] = M[1]/sqrt(M[0]*M[0]+M[1]*M[1]+M[2]*M[2]);
    m[2] = M[2]/sqrt(M[0]*M[0]+M[1]*M[1]+M[2]*M[2]);
    // Update final values for next step

/*********** edition for NDR starts here ***********/
    d_Rmtj = 1/(Gt*(1+(TMR/(TMR+2))*(m[0]*P[0]+m[1]*P[1]+m[2]*P[2]))+Gsi) - R;
/*********** edition for NDR ends here ***********/
    R = 1/(Gt*(1+(TMR/(TMR+2))*(m[0]*P[0]+m[1]*P[1]+m[2]*P[2]))+Gsi);


}
   	if( initial_state ==0){
   		if( R >= Rp*(1+TMR/2)){
   			writeSuccess[this_id*trials_p_thread + i_trial]=1;
   		}
		else{
			writeSuccess[this_id*trials_p_thread + i_trial]=0;
		}
   	}
   	else {
   		if( R <= Rp*(1+TMR/2)){
   			writeSuccess[this_id*trials_p_thread + i_trial]=1;
   		}
		else{
			writeSuccess[this_id*trials_p_thread + i_trial]=0;
		}
   	}

//Recording switching energy and voltage after switching
	g_Energy[this_id * trials_p_thread+i_trial] = energy;
	if(isNDR <=1){
		g_EndVndr[this_id * trials_p_thread+i_trial] = Vndr;
	}
	else{
		g_EndVndr[this_id * trials_p_thread+i_trial] = Vndr+Vmos+V;
	}

}
                  
}


