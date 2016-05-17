//*****************
//This Solver solves the current and voltage of NDR in a series system
// of One NDR, one transistor, and one MTJ. The Bit-line capacitance is
// paralley connected to the point between NDR and transistor.
//#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;
#define N_point 1000 // the number of (V,I) for NDR and transistor (source to drain)

//Read voltage and current characteristics of NDR or Transistor (Vds, Ids)
double* Read_voltage_current(char* file){
    int Npoint = N_point;
    fstream fi;
    fi.open(file,std::fstream::in);
    vector<double> v_v, v_i;
    double readBuffer,voltage,current;
    while(fi >> readBuffer){
        v_v.push_back(readBuffer);
        fi >> readBuffer;
        v_i.push_back(readBuffer);
    }
    fi.close();
    double *VIG = new double[Npoint*3]; //The Voltage and Current and Gradients
//    for(int i=0; i<Npoint; i++){
//        VIG[i] = new double[3];
//    }
    double v_end = v_v[v_v.size() - 1];
    double stepV = v_end/Npoint;
    unsigned int nextI=0; // the index of the voltage in readin data > current dealing V = i*stepV
    cout<<"read in V_I_G."<<endl;
    for (int i = 0; i<Npoint; i++){
        VIG[i*3] = i*stepV;
        unsigned int ind = nextI;
        //Find the next ind with v_v[ind] > i*stepV
        while( ind  <= v_v.size()){
            if(v_v[ind] > i*stepV){
                nextI = ind;
                break;
            }
            ind++;
        }
        //Asign current to V[i][1] and gradient to V[i][2]
        if( ind > 0){
            VIG[i*3+1] = (v_i[ind] - v_i[ind-1])/(v_v[ind] - v_v[ind-1]) * (stepV*i - v_v[ind-1]) + v_i[ind-1];
            VIG[i*3+2] = (v_i[ind] - v_i[ind-1])/(v_v[ind] - v_v[ind-1]) ;
        }
        else{
            VIG[i*3+1] = v_i[ind]/v_v[ind] * stepV*i;
            VIG[i*3+2] = v_i[ind]/v_v[ind];
        }
        if(VIG[i*3+1] < 0)
            VIG[i*3+1] =0;
    //    cout<<VIG[i][0]<<" "<<VIG[i][1]<<" "<<VIG[i][2]<<endl;
    }
    return VIG;
}

//The function returns the current/gradient of NDR or Transistor given voltage, and VI characteristics
// V: voltage, VIG: VIG data base, flag: 1->return current, 2->return gradient
__device__ double IG_V(double V, const double* VIG, int flag){
    double sign = 1;
    if(V < 0){
        sign = -1;
        V = -V;
    }
    //Bisection search
    if( V > VIG[ (N_point-1)*3 ] || V < VIG[0] ){
        //out of bound
        return -32768;
    }
    int begin = 0, end = N_point-1;
    double min_v_step = VIG[3] - VIG[0]; // the precision of voltage
    double diff = VIG[ int( (begin+end)/2)*3] - V;
    while( end > begin+1){
        if(diff > 0){
            end = int((begin+end)/2);
        }
        else{
            begin = int((begin+end)/2);
        }
        diff = VIG[int( (begin+end)/2)*3] - V;
//        cout<<"begin/end:"<<begin<<" / "<<end<<endl;
    }
    if(flag == 1){
        return sign*VIG[int( (begin+end)/2)*3 +  flag];
    }
    else{
        return VIG[int( (begin+end)/2)*3 + flag];
    }

}

//The function returns the voltage given the current
__device__ double V_I(double I, const double* VIG){
    double sign = 1;
    if( I < 0){ // symmetry I-V
        I = -I;
        sign = -1;
    }
    if( I > VIG[ (N_point-1)*3+ 1] ){
        return -32768;
    }
    if( I<VIG[1])
        return 0;

    int begin = 0, end = N_point-1;
    double diff = VIG[int( (begin+end)/2)*3 + 1] - I;
    while( begin < end - 1 ){
        if(diff > 0){
            end = int((begin+end)/2);
        }
        else{
            begin = int((begin+end)/2);
        }
        diff = VIG[int((begin+end)/2)*3 + 1] - I;
    }
    return sign*VIG[int((begin+end)/2)*3 + 0];
}

//This is the stable gillhoff current equation: Vndr + Vt(Indr(Vndr)) + Indr(Vndr) * Rmtj - vdd
__device__ double Equation(double Vndr, const double* VIGndr, const double* VIGmos, double Rmtj, double vdd){
    double I = IG_V(Vndr, VIGndr, 1);
    if( abs(I+32768) <1e-16){
        return -32768;
    }
//    cout<<"I:"<<I<<endl;
//    cout<<"V_mos:"<<V_I(I,VIGmos)<<endl;
    return Vndr + V_I(I, VIGmos) + I*Rmtj - vdd;
}
//Find the peak voltage of NDR
__device__ double Peak_voltage(const double* VIGndr){
    int peak_ind = 0;
    double peak_current = 0;
    for (int i=0; i< N_point; i++){
        if(VIGndr[i*3+1] > peak_current){
            peak_current = VIGndr[i*3+1];
            peak_ind = i;
        }
        if( VIGndr[i*3+1] < 0.5*peak_current){
            break;
        }
    }
    return VIGndr[peak_ind*3];
}

	
//The Solve_stable_vndr solve the stable solution Vndr between (v,i) = (0,0) and (vpeak,Ipeak) from initial condition (0,0), if there is no solution, return -1
__device__ double Solve_stable_vndr(const double* VIGndr, const double* VIGmos, double Rmtj, double vdd){
    //find the peak I and V
    double sign = 1; // sign of output
    if(vdd < 0){
        sign = -1;
        vdd = -vdd;
    }
    double peak_current = 0;
    int peak_ind = 0;
    for (int i=0; i< N_point; i++){
        if(VIGndr[i*3+1] > peak_current){
            peak_current = VIGndr[i*3+1];
            peak_ind = i;
        }
        if( VIGndr[i*3+1] < 0.5*peak_current){
            break;
        }
    }
    //cout<<"peak voltage is:"<<sign*VIGndr[peak_ind][0]<<endl;
//Using bisection to solve the equation Vndr + Vt(Indr(Vndr)) + Indr(Vndr) * Rmtj - vdd =0
    double v_begin = 0, v_end = (VIGndr[peak_ind*3+0] < vdd) ? VIGndr[peak_ind*3+0] : vdd;
    double e_begin =  Equation(v_begin,VIGndr,VIGmos,Rmtj,vdd);
    double e_end = Equation(v_end,VIGndr,VIGmos,Rmtj,vdd);
    if( e_begin > 0 || fabs(e_begin+32768) < 1e-16 || fabs(e_end+32768)<1e-16 ) {
        // No solution
        return -32768;
    }
    if( e_end < 0){
        if(vdd > VIGndr[peak_ind*3]){
            v_begin = VIGndr[peak_ind*3+0];
            v_end = (VIGndr[(N_point-1)*3 +0] < vdd) ? VIGndr[(N_point-1)*3+0] : vdd;
        }
        else{
            return -32768;
        }
    }
    double stepV = VIGndr[3] - VIGndr[0];
    //cout<<"Hellow world1"<<endl;
    double e_middle = Equation((v_end+v_begin)/2,VIGndr,VIGmos,Rmtj,vdd); 
    while(v_begin < v_end - stepV){
        //cout<<"v_begin/v_end: "<<v_begin<<" "<<v_end<<endl;
        //cout<<"e_middle:"<<e_middle<<endl;
        if(e_middle < 0 ){
            v_begin = (v_begin + v_end)/2;
        }
        else{
            v_end = (v_begin + v_end)/2;
        }
        e_middle = Equation((v_end+v_begin)/2,VIGndr,VIGmos,Rmtj,vdd);
    }
    return sign*(v_begin + v_end)/2;
}
