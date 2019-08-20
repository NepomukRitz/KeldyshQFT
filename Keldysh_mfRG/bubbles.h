#pragma clang diagnostic push
#pragma ide diagnostic ignored "bugprone-narrowing-conversions"

//
// Created by E.Walter on 8/1/19.
//

#ifndef KELDYSH_MFRG_BUBBLES_H
#define KELDYSH_MFRG_BUBBLES_H

#include <algorithm>
#include <gsl/gsl_integration.h>

#include "selfenergy.h"
#include "propagator.h"
#include "data_structures.h"
#include "frequency_grid.h"

using namespace std;


/******************************************** BUBBLE INTEGRATION FUNCTIONS ********************************************/

//a-bubble:
template<typename  T1,typename  T2>
struct abubble_params{
    int red_side;
    int map1;
    int map2;
    double Lambda;
    T1& vert1;
    T2& vert2;
    char ptype1;
    char ptype2;
    SelfEnergy<comp>& selfenergy;
    SelfEnergy<comp>& diffselfenergy;

    double w;double v1; double v2;
    char h;//class of bubble (K = K1, L=K2, M = K2b, R = K3)

};

template<typename T1,typename T2>
double abubble_re(double u, void * p){
    auto * params
            = static_cast< struct abubble_params<T1,T2> *>(p);
    char red_side = (params ->red_side);
    char map1 = (params ->map1);
    char map2 = (params ->map2);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
    SelfEnergy<comp>& selfenergy = (params->selfenergy);
    SelfEnergy<comp>& diffselfenergy = (params->diffselfenergy);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double w = (params->w);
    double v1 = (params->v1);
    double v2 = (params->v2);
    char h = (params->h);

    double val = real(vert1.vvalsmooth(red_side,map1,a,b,c,w,u,v2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,w,v1,u,'u',2,h) *  propag(Lambda,u-w/2,selfenergy,diffselfenergy,ptype1) * propag(Lambda,u+w/2,selfenergy,diffselfenergy,ptype2) );
    return (1./(2*pi)*val);
}

template<typename T1,typename T2>
double agreensfunc_re(double u, void * p){
    auto * params
            = static_cast< struct abubble_params<T1,T2> *>(p);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
    SelfEnergy<comp>& selfenergy = (params->selfenergy);
    SelfEnergy<comp>& diffselfenergy = (params->diffselfenergy);
    double Lambda = (params->Lambda);

    double  w = (params->w);


    double val = real( propag(Lambda,u-w/2,selfenergy,diffselfenergy,ptype1) * propag(Lambda,u+w/2,selfenergy,diffselfenergy,ptype2) );
    return (1./(2*pi)*val);
}

template<typename T1,typename T2>
double abubble(int red_side,int map1,int map2,gsl_integration_workspace* w,double Lambda, T1& vert1,int a, int b, int c, T2& vert2,int d, int e, int f,char p1, char p2, SelfEnergy<comp>& se, SelfEnergy<comp>& dse, double u,double w1, double w2, char h){//bubble in s-channel: specified by s: bos. frequency and two external ferm freqs.

    double abs_error=1e-2, rel_error=1e-4;

    double abs_error_bare=1e-4,rel_error_bare=1e-4;
    double B=0.;

    if((REG ==1 && p1 =='s')){//if first propagator is single scale propagator, proportional to delta peak: no integration
        B = real(1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda+u/2,'u',2,h)* propag(Lambda,-Lambda,se,dse,p1) * propag(Lambda,-Lambda+u,se,dse,p2)) ;
        B += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda+u/2,'u',2,h) *  propag(Lambda,Lambda,se,dse,p1) * propag(Lambda,Lambda+u,se,dse,p2)) ;
    }
    else if((REG ==1 && p2 =='s') ){//if second propagator is single scale propagator, proportional to delta peak: no integration

        B =  real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda-u/2,'u',2,h) *  propag(Lambda,-Lambda-u,se,dse,p1) * propag(Lambda,-Lambda,se,dse,p2)) ;
        B += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda-u/2,'u',2,h) *  propag(Lambda,Lambda-u,se,dse,p1) * propag(Lambda,Lambda,se,dse,p2)) ;
    }





    else{
        double limit_low,limit_low1,limit_low2;
        double limit_up,limit_up1,limit_up2;



        double vert1_const,vert2_const;




        int mode =1;

        if(h=='R'){


            //upper bound:

            //conditions from vertex 1 (primed):
            double v1_K1_t_up = bfreqs[nw1-1]+w2;
            double v1_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+u},compare);
            double v1_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-u},compare);
            double v1_K3_t_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2+u,2*ffreqs[(nw+nw3)/2-1]-w2-u},compare);
            double v1_K1_s_up= bfreqs[nw1-1]-w2;
            double v1_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+u},compare);
            double v1_K2b_s_up =  min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-u},compare);
            double v1_K3_s_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2-u,2*ffreqs[(nw+nw3)/2-1]+w2+u},compare);

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_up = bfreqs[nw1-1]+w1;
            double v2_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1+u},compare);
            double v2_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1-u},compare);
            double v2_K3_t_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw1+nw3)/2-1]-w1-u,2*ffreqs[(nw1+nw3)/2-1]-w1+u},compare);
            double v2_K1_s_up = bfreqs[nw1-1]-w1;
            double v2_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1-u},compare);
            double v2_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1+u},compare);
            double v2_K3_s_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw+nw3)/2-1]+w1+u,2*ffreqs[(nw+nw3)/2-1]+w1-u},compare);
            //conditions that is equal for both vertices (u-channel contributions). Note that K_{a,1} is constant
            double v_K2_u_up = ffreqs[(nw1+nw2)/2-1];
            double v_K3_u_up = ffreqs[(nw1+nw3)/2-1];



            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_t_low = bfreqs[0]+w2;
            double v1_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+u},compare);
            double v1_K2b_t_low =max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-u},compare);
            double v1_K3_t_low = max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2+u,2*ffreqs[(nw-nw3)/2]-w2-u},compare);
            double v1_K1_s_low = bfreqs[0]-w2;
            double v1_K2_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+u},compare);
            double v1_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-u},compare);
            double v1_K3_s_low = max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2-u,2*ffreqs[(nw-nw3)/2]+w2+u},compare);


            //conditions from vertex 2 (unprimed):
            double v2_K1_t_low = bfreqs[0]+w1;
            double v2_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1+u},compare);
            double v2_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1-u},compare);
            double v2_K3_t_low =  max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw1-nw3)/2]-w1-u,2*ffreqs[(nw1-nw3)/2]-w1+u},compare);
            double v2_K1_s_low = bfreqs[0]-w1;
            double v2_K2_s_low =  max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1-u},compare);
            double v2_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1+u},compare);
            double v2_K3_s_low =  max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw-nw3)/2]+w1+u,2*ffreqs[(nw-nw3)/2]+w1-u},compare);
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_u_low = ffreqs[(nw1-nw2)/2];
            double v_K3_u_low = ffreqs[(nw1-nw3)/2];


            vector<double> upperK3_v1{v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_K3_s_up,v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_K3_t_up,v_K2_u_up,v_K3_u_up};
            sort(upperK3_v1.begin(),upperK3_v1.end());
            vector<double> lowerK3_v1{v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_K3_s_low,v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_K3_t_low,v_K2_u_low,v_K3_u_low};
            sort(lowerK3_v1.begin(),lowerK3_v1.end());

            vector<double> upperK3_v2{v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_K3_s_up,v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_K3_t_up,v_K2_u_up,v_K3_u_up};
            sort(upperK3_v2.begin(),upperK3_v2.end());
            vector<double> lowerK3_v2{v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_K3_s_low,v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_K3_t_low,v_K2_u_low,v_K3_u_low};
            sort(lowerK3_v2.begin(),lowerK3_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,u,wlimit,w2,'u',1,h);
            vert2_const = vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,wlimit,'u',2,h);

            for(int i=upperK3_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK3_v2[i],'u',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK3_v2[i],w2,'u',1,h))>0 ){limit_up2 = upperK3_v2[i];break;}
                if(i==0){mode2_up=0;}
            }

            for(int i=upperK3_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK3_v1[i],w2,'u',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK3_v1[i],'u',2,h))>0 ){limit_up1 = upperK3_v1[i];break;}
                if(i==0){mode1_up=0;}
            }

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;}



            for(int i=0; i< lowerK3_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK3_v1[i],w2,'u',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK3_v1[i],'u',2,h))>0){limit_low1 = lowerK3_v1[i];break;}
                if(i==lowerK3_v1.size()-1){mode1_low=0;}
            }

            for(int i=0; i< lowerK3_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK3_v2[i],'u',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK3_v2[i],w2,'u',1,h))>0){limit_low2 = lowerK3_v2[i];break;}
                if(i==lowerK3_v2.size()-1){mode2_low=0;}
            }


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;}
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;}}


        else if(h=='L'){ //constant vertices



            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_u_low = ffreqs[(nw1-nw2)/2];
            double v_K3_u_low = ffreqs[(nw1-nw3)/2];
            //conditions that is equal for both vertices (u-channel contributions). Note that K_{a,1} is constant
            double v_K2_u_up = ffreqs[(nw1+nw2)/2-1];
            double v_K3_u_up = ffreqs[(nw1+nw3)/2-1];

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_low = bfreqs[0]+w1;
            double v2_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1+u},compare);
            double v2_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1-u},compare);
            double v2_K3_t_low =  max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw1-nw3)/2]-w1-u,2*ffreqs[(nw1-nw3)/2]-w1+u},compare);
            double v2_K1_s_low = bfreqs[0]-w1;
            double v2_K2_s_low =  max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1-u},compare);
            double v2_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1+u},compare);
            double v2_K3_s_low =  max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw-nw3)/2]+w1+u,2*ffreqs[(nw-nw3)/2]+w1-u},compare);

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_up = bfreqs[nw1-1]+w1;
            double v2_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1+u},compare);
            double v2_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1-u},compare);
            double v2_K3_t_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw1+nw3)/2-1]-w1-u,2*ffreqs[(nw1+nw3)/2-1]-w1+u},compare);
            double v2_K1_s_up = bfreqs[nw1-1]-w1;
            double v2_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1-u},compare);
            double v2_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1+u},compare);
            double v2_K3_s_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw+nw3)/2-1]+w1+u,2*ffreqs[(nw+nw3)/2-1]+w1-u},compare);

            vector<double> upperK2_v1{v_K2_u_up};
            sort(upperK2_v1.begin(),upperK2_v1.end());
            vector<double> lowerK2_v1{v_K2_u_low};
            sort(lowerK2_v1.begin(),lowerK2_v1.end());

            vector<double> upperK2_v2{v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_K3_s_up,v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_K3_t_up,v_K2_u_up,v_K3_u_up};
            sort(upperK2_v2.begin(),upperK2_v2.end());
            vector<double> lowerK2_v2{v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_K3_s_low,v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_K3_t_low,v_K2_u_low,v_K3_u_low};
            sort(lowerK2_v2.begin(),lowerK2_v2.end());


            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,u,wlimit,wlimit,'u',1,h);
            vert2_const = vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,wlimit,'u',2,h);


            for(int i=upperK2_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK2_v2[i],'u',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK2_v2[i],w2,'u',1,h))>0 ){limit_up2 = upperK2_v2[i];break;}
                if(i==0){mode2_up=0;}

            }

            for(int i=upperK2_v1.size()-1; i>-1 ; i--){

                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK2_v1[i],w2,'u',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK2_v1[i],'u',2,h))>0 ){limit_up1 = upperK2_v1[i];break;}
                if(i==0){mode1_up=0;}
            }

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;}

            for(int i=0; i< lowerK2_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK2_v1[i],w2,'u',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK2_v1[i],'u',2,h))>0){limit_low1 = lowerK2_v1[i];break;}
                if(i==lowerK2_v1.size()-1){mode1_low=0;}
            }

            for(int i=0; i< lowerK2_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK2_v2[i],'u',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK2_v2[i],w2,'u',1,h))>0){limit_low2 = lowerK2_v2[i];break;}
                if(i==lowerK2_v2.size()-1){mode2_low=0;}
            }


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;}

            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;}
        }

        else if(h=='M'){//in this case, only Gamma' (vert1) sets the limits

            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_u_low = ffreqs[(nw1-nw2)/2];
            double v_K3_u_low = ffreqs[(nw1-nw3)/2];

            //conditions that is equal for both vertices (u-channel contributions). Note that K_{a,1} is constant
            double v_K2_u_up = ffreqs[(nw1+nw2)/2-1];
            double v_K3_u_up = ffreqs[(nw1+nw3)/2-1];


            //conditions from vertex 1 (primed):
            double v1_K1_t_up = bfreqs[nw1-1]+w2;
            double v1_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+u},compare);
            double v1_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-u},compare);
            double v1_K3_t_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2+u,2*ffreqs[(nw+nw3)/2-1]-w2-u},compare);
            double v1_K1_s_up= bfreqs[nw1-1]-w2;
            double v1_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+u},compare);
            double v1_K2b_s_up =  min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-u},compare);
            double v1_K3_s_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2-u,2*ffreqs[(nw+nw3)/2-1]+w2+u},compare);

            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_t_low = bfreqs[0]+w2;
            double v1_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+u},compare);
            double v1_K2b_t_low =max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-u},compare);
            double v1_K3_t_low = max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2+u,2*ffreqs[(nw-nw3)/2]-w2-u},compare);
            double v1_K1_s_low = bfreqs[0]-w2;
            double v1_K2_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+u},compare);
            double v1_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-u},compare);
            double v1_K3_s_low = max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2-u,2*ffreqs[(nw-nw3)/2]+w2+u},compare);


            vector<double> upperK2b_v1{v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_K3_s_up,v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_K3_t_up,v_K2_u_up,v_K3_u_up};
            sort(upperK2b_v1.begin(),upperK2b_v1.end());
            vector<double> lowerK2b_v1{v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_K3_s_low,v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_K3_t_low,v_K2_u_low,v_K3_u_low};
            sort(lowerK2b_v1.begin(),lowerK2b_v1.end());

            vector<double> upperK2b_v2{v_K2_u_up};
            sort(upperK2b_v2.begin(),upperK2b_v2.end());
            vector<double> lowerK2b_v2{v_K2_u_low};
            sort(lowerK2b_v2.begin(),lowerK2b_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;

            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,u,wlimit,w2,'u',1,h);
            vert2_const =vert2.vvalsmooth(red_side,map2,d,e,f,u,wlimit,wlimit,'u',2,h);

            for(int i=upperK2b_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK2b_v2[i],'u',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK2b_v2[i],w2,'u',1,h))>0 ){limit_up2 = upperK2b_v2[i];break;}
                if(i==0){mode2_up=0;}
            }

            for(int i=upperK2b_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK2b_v1[i],w2,'u',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK2b_v1[i],'u',2,h))>0 ){limit_up1 = upperK2b_v1[i];break;}
                if(i==0){mode1_up=0;}
            }


            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;}


            for(int i=0; i< lowerK2b_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK2b_v1[i],w2,'u',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK2b_v1[i],'u',2,h))>0){limit_low1 = lowerK2b_v1[i];break;}
                if(i==lowerK2b_v1.size()-1){mode1_low=0;}
            }

            for(int i=0; i< lowerK2b_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK2b_v2[i],'u',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK2b_v2[i],w2,'u',1,h))>0){limit_low2 = lowerK2b_v2[i];break;}
                if(i==lowerK2b_v2.size()-1){mode2_low=0;}
            }


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;}}

        else if(h=='K'){ //in this case the only contribution to the non-constat part comes from K2,u


            double v_K2_u_low = ffreqs[(nw1-nw2)/2];
            double v_K2_u_up = ffreqs[(nw1+nw2)/2-1];

            vector<double> upperK1_v1{v_K2_u_up};
            sort(upperK1_v1.begin(),upperK1_v1.end());
            vector<double> lowerK1_v1{v_K2_u_low};
            sort(lowerK1_v1.begin(),lowerK1_v1.end());

            vector<double> upperK1_v2{v_K2_u_up};
            sort(upperK1_v2.begin(),upperK1_v2.end());
            vector<double> lowerK1_v2{v_K2_u_low};
            sort(lowerK1_v2.begin(),lowerK1_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,u,wlimit,wlimit,'u',1,h);
            vert2_const = vert2.vvalsmooth(red_side,map2,d,e,f,u,wlimit,wlimit,'u',2,h);



            for(int i=upperK1_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK1_v2[i],'u',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK1_v2[i],w2,'u',1,h))>0 ){limit_up2 = upperK1_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK1_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,upperK1_v1[i],w2,'u',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,upperK1_v1[i],'u',2,h))>0 ){limit_up1 = upperK1_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            for(int i=0; i< lowerK1_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK1_v1[i],w2,'u',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK1_v1[i],'u',2,h))>0){limit_low1 = lowerK1_v1[i];break;};
                if(i==lowerK1_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK1_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,lowerK1_v2[i],'u',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,u,lowerK1_v2[i],w2,'u',1,h))>0){limit_low2 = lowerK1_v2[i];break;};
                if(i==lowerK1_v2.size()-1){mode2_low=0;};
            };

            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};};

        double result_re=0, error_re=0;
        double resultgfl_re=0, errorgfl_re=0;
        double resultgfu_re=0, errorgfu_re=0;

        struct abubble_params<T1,T2> params = {red_side,map1,map2,Lambda,vert1, vert2,p1,p2,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};
        gsl_function F;

        if(mode==0 && abs(vert1_const * vert2_const)>0){

            F.function = &agreensfunc_re<T1,T2>;

            if(REG==1 && p1=='g' && p2 =='g'){
                F.params = &params;


                double upper = abs(u/2) + Lambda;

                //integrand is even: sufficient to integrate positive part and multiply result by 2



                gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                      w, &resultgfu_re, &errorgfu_re);


                double compl_real = 2*resultgfu_re *vert1_const * vert2_const;

                B= compl_real ;
            }


            else if(REG==1 && ((p1=='k' &&p2=='g') ||(p1=='g'&& p2 =='k'))){
                double lower = -abs(u/2)-Lambda;
                double upper = abs(u/2) + Lambda;

                char p1_new,p2_new;
                double dse_cutoff_up=0, dse_cutoff_low=0;//cutoff due to finite number of saved value in differentiated self energy
                double singsc=0;




                if(p1=='k'){
                    singsc =  real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda+u/2,w2,'u',1,h)* vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda+u/2,'u',2,h)* propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+u,se,dse,p2)) ;
                    singsc += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u, Lambda+u/2,w2,'u',1,h)* vert2.vvalsmooth(red_side,map2,d,e,f,u,w1, Lambda+u/2,'u',2,h)* propag(Lambda, Lambda,se,dse,'s') * propag(Lambda, Lambda+u,se,dse,p2)) ;
                    p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]+u/2;dse_cutoff_low= ffreqs[0]+u/2;}

                else if(p2=='k'){
                    singsc =  real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda-u/2,'u',2,h) *  propag(Lambda,-Lambda-u,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s') );
                    singsc += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u, Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1, Lambda-u/2,'u',2,h) *  propag(Lambda, Lambda-u,se,dse,p1) * propag(Lambda, Lambda,se,dse,'s')) ;
                    p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]-u/2;dse_cutoff_low= ffreqs[0]-u/2;};

                struct abubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(dse_cutoff_low <lower && dse_cutoff_up <=upper){//if highest saved dse-value if lower than upper heavyside-cutoff

                    if(dse_cutoff_up>lower){dse_cutoff_up=lower;}

                    F.params = &params_mod;

                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                    double compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    B= singsc+ compl_real ;
                }

                else if(dse_cutoff_low >=lower && dse_cutoff_up > upper){//if lowest saved dse-value if higher than lower heavyside-cutoff

                    if(dse_cutoff_low<upper){dse_cutoff_low=upper;};

                    F.params = &params_mod;

                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                    double compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    B=singsc + compl_real ;
                }
                else if (dse_cutoff_low<lower && dse_cutoff_up > upper){
                    F.params = &params_mod;


                    gsl_integration_qag(&F,dse_cutoff_low,lower,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    gsl_integration_qag(&F,upper,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfu_re, &errorgfu_re);

                    double compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                    B=singsc + compl_real ;
                }}


            else if(REG==2 && (p1=='g' &&p2=='g')){

                F.params = &params;


                gsl_integration_qagiu(&F,0,abs_error_bare,rel_error_bare,1500,
                                      w, &result_re, &errorgfu_re);


                double compl_real = 2*result_re *vert1_const * vert2_const;

                B=compl_real;


            }
            else if(REG==2 && ((p1=='g' &&p2=='s') ||(p1=='s'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially
                double bound_low=0, bound_up=0;
                if(p1=='s'){bound_low = -7*Lambda+u/2; bound_up = 7*Lambda+u/2;}
                else if(p2=='s'){bound_low = -7*Lambda-u/2; bound_up = 7*Lambda-u/2;}
                F.params = &params;
                gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);

                double compl_real =(result_re) *vert1_const * vert2_const;

                B=compl_real;

            }

            else if(REG==2 && ((p1=='g' &&p2=='k') ||(p1=='k'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially

                double bound_low=0, bound_up=0, dse_cutoff_low=0, dse_cutoff_up=0;
                char p1_new, p2_new, p1_singsc,p2_singsc;
                if(p1=='k'){
                    bound_low = -7*Lambda+u/2; bound_up = 7*Lambda+u/2;
                    dse_cutoff_low = ffreqs[0]+u/2; dse_cutoff_up = ffreqs[nw-1]+u/2;
                    p1_new = 'e';
                    p1_singsc = 's';
                    p2_new = 'g';
                    p2_singsc = 'g';
                }
                else if(p2=='k'){
                    bound_low = -7*Lambda-u/2; bound_up = 7*Lambda-u/2;
                    dse_cutoff_low = ffreqs[0]-u/2; dse_cutoff_up = ffreqs[nw-1] -u/2;
                    p1_new = 'g';
                    p1_singsc = 'g';
                    p2_new = 'e';
                    p2_singsc = 's';
                }

                //first compute single scale contribution
                struct abubble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_singsc,p2_singsc,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};
                F.params = &params_singsc;
                gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);

                //add katanin extension
                struct abubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                F.params = &params_mod;//add integration of Katanin extension
                double result_re2=0,error_re2=0;
                gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error, rel_error,1500,2,
                                    w, &result_re2, &error_re2);

                double compl_real =(result_re + result_re2 ) *vert1_const * vert2_const;
                B=compl_real;
            }
        }

        else if(mode==1){


            if(REG==1 && p1=='g' && p2=='g'){
                double lhs=0, rhs=0;//result on left and right side of sharp cutoff
                double lower = -abs(u/2)-Lambda;
                double upper = abs(u/2) + Lambda;
                if(lower <= limit_low){//left side of sharp cutoff

                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,lower,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);

                        double compl_real =(resultgfl_re  )*vert1_const * vert2_const;

                        lhs=compl_real;

                    }}
                else if(lower <= limit_up){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);
                        compl_real =(resultgfl_re  )*vert1_const * vert2_const;

                    }

                    F.params = &params;
                    F.function = &abubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,lower,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    lhs=result_re + compl_real;

                }

                else if(lower > limit_up){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,lower,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };
                    F.params = &params;
                    F.function = &abubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    lhs=result_re + compl_real;


                };

                if(upper >= limit_up){//right side of sharp cutoff

                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                        double compl_real = (resultgfu_re  )*vert1_const * vert2_const;

                        rhs=compl_real;

                    };}
                else if(upper >= limit_low){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                        compl_real =(resultgfu_re  )*vert1_const * vert2_const;
                    };
                    F.params = &params;
                    F.function = &abubble_re<T1,T2>;

                    gsl_integration_qag(&F,upper,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    rhs=result_re + compl_real;

                }

                else if (upper < limit_low){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,upper,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);

                        compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                    };
                    F.params = &params;
                    F.function = &abubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    rhs=result_re + compl_real;


                };

                B= rhs + lhs;


            }
            else if(REG==1 && ((p1=='k' &&p2=='g') || (p1=='g'&& p2 =='k'))){

                double lower = -abs(u/2)-Lambda;
                double upper = abs(u/2) + Lambda;

                char p1_new,p2_new;
                double dse_cutoff_up=0, dse_cutoff_low=0;//cutoff due to finite number of saved value in differentiated self energy
                double singsc=0;

                if(p1=='k'){
                    singsc = real( 1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda+u/2,'u',2,h) *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+u,se,dse,p2) );
                    singsc += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda+u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,Lambda+u/2,'u',2,h) *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,Lambda+u,se,dse,p2)) ;
                    p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]+u/2;dse_cutoff_low= ffreqs[0]+u/2;}

                else if(p2=='k'){
                    singsc =  real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,-Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,u,w1,-Lambda-u/2,'u',2,h) *  propag(Lambda,-Lambda-u,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
                    singsc += real(1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,u,Lambda-u/2,w2,'u',1,h) * vert2.vvalsmooth(d,e,f,u,w1,Lambda-u/2,'u',2,h) *  propag(Lambda,Lambda-u,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
                    p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]-u/2;dse_cutoff_low= ffreqs[0]-u/2;};
                double lhs=0, rhs=0;//result on left and right side of sharp cutoff

                struct abubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(lower <= limit_low && abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){//left side of sharp cutoff

                    double up_eff=0;
                    if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                    F.params = &params_mod;
                    F.function = &agreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low ,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                    double compl_real = (resultgfl_re )*vert1_const * vert2_const;
                    lhs=compl_real;
                }
                else if(lower <= limit_up ){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    };
                    if(dse_cutoff_low <lower){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                        F.params = &params_mod;
                        F.function = &abubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                    };

                    lhs=(result_re + compl_real);

                }

                else if(lower > limit_up){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    };


                    if(dse_cutoff_low <limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &abubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                    };


                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                        double up_eff=0;
                        if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                        F.params = &params_mod;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real += (resultgfl_re )*vert1_const * vert2_const;
                    };

                    lhs=(result_re + compl_real);


                };

                if(upper >= limit_up && abs(vert1_const * vert2_const)>0 && dse_cutoff_up >upper){//right side of sharp cutoff

                    double low_eff=0;
                    if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                    F.params = &params_mod;
                    F.function = &agreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfu_re, &errorgfu_re);
                    double compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    rhs=compl_real;
                }
                else if(upper >= limit_low ){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up >limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfu_re  )*vert1_const * vert2_const;

                    };
                    if(dse_cutoff_up > upper){
                        double low_eff=0;
                        if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &abubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                    };

                    rhs=(result_re + compl_real);

                }

                else if(upper < limit_low){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > upper){
                        double low_eff=0;
                        if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;

                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real =(resultgfu_re  )*vert1_const * vert2_const;

                    };


                    if(dse_cutoff_up >limit_low){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &abubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                    };


                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};

                        F.params = &params_mod;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real +=(resultgfu_re )*vert1_const * vert2_const;

                    };

                    rhs=(result_re + compl_real);


                };

                B=rhs +lhs+ singsc;
            }

            else if(REG==2 && p1=='g' && p2=='g'){

                F.function = &abubble_re<T1,T2>;
                F.params = &params;

                gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);

                //integration of greens function outside of w-dependent interval

                double compl_real=0;
                if(abs(vert1_const * vert2_const)>0){
                    F.function = agreensfunc_re<T1,T2>;
                    gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfl_re, &errorgfl_re);

                    F.function = &agreensfunc_re<T1,T2>;
                    gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfu_re, &errorgfu_re);

                    compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                };
                B=(result_re + compl_real );
            }

            else if(REG==2 && ((p1=='g' && p2=='s') || (p1=='s' && p2=='g'))){


                double bound_low=0, bound_up=0;

                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;

                if(p1=='s'){bound_low = -7*Lambda+u/2; bound_up = 7*Lambda+u/2;}
                else if(p2=='s'){bound_low = -7*Lambda-u/2; bound_up = 7*Lambda-u/2;}

                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };




                    B=result_re + compl_real ;
                }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real =(resultgfl_re )*vert1_const * vert2_const;

                    };


                    B=(result_re + compl_real );
                }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){

                        F.params = &params;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);


                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };



                    B=result_re + compl_real ;
                }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    B=( compl_real );
                }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    B=(result_re  );
                }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };

                    B=( compl_real );
                }

            }

            else if(REG==2 && ((p1=='g' && p2=='k') || (p1=='k' && p2=='g'))){

                double singsc=0;
                double bound_low=0, bound_up=0;
                //FIRST COMPUTE THE SINGLE SCALE CONTRIBUTION
                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;
                char p1_new,p2_new;
                if(p1=='k'){
                    bound_low = -7*Lambda+u/2; bound_up = 7*Lambda+u/2;
                    p1_new = 's'; p2_new = 'g';}
                else if(p2=='k'){bound_low = -7*Lambda-u/2; bound_up = 7*Lambda-u/2;
                    p1_new = 'g'; p2_new = 's';
                }


                struct abubble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };




                    singsc=(result_re + compl_real );
                }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    singsc=result_re + compl_real ;
                }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){

                        F.params = &params_singsc;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };



                    singsc=(result_re + compl_real);
                }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    singsc=( compl_real );
                }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    singsc=(result_re  );
                }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    };

                    singsc=( compl_real );
                }
                //ADD THE CONTRIBUTION FROM THE KATANON EXTENSION:
                double kat=0;

                //reset m-values
                m1=0;//turns on and off the integration in the three parts
                m2=0;
                m3=0;


                if(p1=='k'){
                    bound_low = ffreqs[0]+u/2; bound_up = ffreqs[nw-1]+u/2;//the bounds are the boundaries of the dse-grid
                    p1_new = 'e';
                    p2_new = 'g';}
                else if(p2=='k'){
                    bound_low = ffreqs[0]-u/2; bound_up = ffreqs[nw-1]-u/2;
                    p1_new = 'g';
                    p2_new = 'e';}

                struct abubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  u, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfl_re +resultgfu_re)*vert1_const * vert2_const;

                    };




                    kat=(result_re + compl_real );
                }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    kat=(result_re + compl_real );
                }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){


                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };



                    kat=(result_re + compl_real );
                }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);


                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    kat=( compl_real );
                }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &abubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    kat=result_re  ;
                }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                    //   double compl_imag;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &agreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    };

                    kat=( compl_real );
                }
                B=(kat + singsc);

            };


        };



        //cout << " ok" << endl;



    };
    if(abs(B)<1e-20){B=0;};

    return B;

}



//p-bubble:
template<typename  T1,typename  T2>
struct pbubble_params{


    int red_side;
    int map1;
    int map2;
    double Lambda;
    T1& vert1;
    T2& vert2;
    char ptype1;
    char ptype2;
    SelfEnergy<comp>& selfenergy;
    SelfEnergy<comp>& diffselfenergy;


    int a; int b; int c;

    int d; int e; int f;

    double s;double w1; double w2;
    char h;

};

template<typename T1,typename T2>
double pbubble_re(double w, void * p){
    struct pbubble_params<T1,T2> * params
            = static_cast< struct pbubble_params<T1,T2> *>(p);

    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
    SelfEnergy<comp>& selfenergy = (params->selfenergy);
    SelfEnergy<comp>& diffselfenergy = (params->diffselfenergy);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);

    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double s = (params->s);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params->h);
    double val = real(1./2 * vert1.vvalsmooth(red_side,map1,a,b,c,s,w,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,w,'s',2,h)  *  propag(Lambda,w+s/2,selfenergy, diffselfenergy,ptype1) * propag(Lambda,s/2-w,selfenergy,diffselfenergy,ptype2) );
    return (1./(2*pi)*val);
}


template<typename T1,typename T2>
double pgreensfunc_re(double w, void * p){//intergrates only propagators - used if vertices are constant in respective integration interval
    struct pbubble_params<T1,T2> * params
            = static_cast< struct pbubble_params<T1,T2> *>(p);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
    SelfEnergy<comp>& selfenergy = (params->selfenergy);
    SelfEnergy<comp>& diffselfenergy = (params->diffselfenergy);
    double Lambda = (params->Lambda);

    double s = (params->s);

    double val = real(1./2 * propag(Lambda,w+s/2,selfenergy, diffselfenergy,ptype1) * propag(Lambda,s/2-w,selfenergy,diffselfenergy,ptype2) );

    return (1./(2*pi)*val);
}


template<typename T1,typename T2>
double pbubble(int red_side,int map1, int map2, gsl_integration_workspace* w, double Lambda, T1& vert1,int a, int b, int c, T2& vert2,int d, int e, int f,char p1, char p2, SelfEnergy<comp>& se, SelfEnergy<comp>& dse, double s,double w1, double w2, char h){//bubble in s-channel: specified by s: bos. frequency and two external ferm freqs.
    double abs_error=1e-2, rel_error=1e-4;
    double abs_error_bare=1e-4,rel_error_bare=1e-4;
    //Note: factor 1/2 is included due to indistiguishability of propagators in s-bubble

    double B =0;

    if((REG ==1 && p1 =='s') ){//if first propagator is single scale propagator, proportional to delta peak: no integration
        B =  real(1./2 * 1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,-Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,-Lambda-s/2,'s',2,h)  *  propag(Lambda,-Lambda,se,dse,p1) * propag(Lambda,s+Lambda,se,dse,p2)) ;
        B += real(1./2 *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,s,Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,Lambda-s/2,'s',2,h)   *  propag(Lambda,Lambda,se,dse,p1) * propag(Lambda,s-Lambda, se,dse, p2));
    }
    else if((REG ==1 && p2 =='s') ){//if second propagator is single scale propagator, proportional to delta peak: no integration
        B =   real(1./2 * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2+Lambda,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2+Lambda,'s',2,h)  *  propag(Lambda,s+Lambda,se,dse,p1) * propag(Lambda,-Lambda,se,dse,p2)) ;
        B += real(1./2 *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2-Lambda,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2-Lambda,'s',2,h)  *  propag(Lambda,s-Lambda,se,dse,p1) * propag(Lambda,Lambda,se,dse,p2)) ;

    }


    else{

        double limit_low,limit_low1,limit_low2;
        double limit_up,limit_up1,limit_up2;



        double vert1_const,vert2_const;



        int mode =1;

        if(h=='R'){
            //upper bound:

            //conditions from vertex 1 (primed):
            double v1_K1_t_up = bfreqs[nw1-1]+w2;
            double v1_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-s},compare);//choose the minimum since beyond this value, this vertex cannot contribute any more
            double v1_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+s},compare);
            double v1_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2-s,2*ffreqs[(nw+nw3)/2-1]-w2+s},compare);
            double v1_K1_u_up= bfreqs[nw1-1]-w2;
            double v1_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-s},compare);
            double v1_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+s},compare);
            double v1_R_u_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2-s,2*ffreqs[(nw+nw3)/2-1]+w2+s},compare);

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_up = bfreqs[nw1-1]+w1;
            double v2_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw+nw2)/2-1]-w1-s},compare);
            double v2_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw+nw2)/2-1]-w1+s},compare);
            double v2_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw+nw3)/2-1]-w1+s,2*ffreqs[(nw+nw3)/2-1]-w1-s},compare);
            double v2_K1_u_up = bfreqs[nw1-1]-w1;
            double v2_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1+s},compare);
            double v2_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1-s},compare);
            double v2_R_u_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw+nw3)/2-1]+w1+s,2*ffreqs[(nw+nw3)/2-1]+w1-s},compare);
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_s_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_s_up = ffreqs[(nw1+nw3)/2-1];



            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_t_low = bfreqs[0]+w2;
            double v1_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-s},compare);//choose the minimum since bezond this value, this vertex cannot contribute any more
            double v1_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+s},compare);
            double v1_R_t_low = max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2-s,2*ffreqs[(nw-nw3)/2]-w2+s},compare);
            double v1_K1_u_low = bfreqs[0]-w2;
            double v1_K2_u_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-s},compare);
            double v1_K2b_u_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+s},compare);
            double v1_R_u_low = max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2-s,2*ffreqs[(nw-nw3)/2]+w2+s},compare);

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_low =  bfreqs[0]+w1;
            double v2_K2_t_low = max({bfreqs[(nw-+nw2)/2]+w1,2*ffreqs[(nw-nw2)/2]-w1-s},compare);
            double v2_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw-nw2)/2]-w1+s},compare);
            double v2_R_t_low = max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw-nw3)/2]-w1+s,2*ffreqs[(nw-nw3)/2]-w1-s},compare);
            double v2_K1_u_low = bfreqs[0]-w1;
            double v2_K2_u_low =  max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1+s},compare);
            double v2_K2b_u_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1-s},compare);
            double v2_R_u_low = max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw-nw3)/2]+w1+s,2*ffreqs[(nw-nw3)/2]+w1-s},compare);
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_s_low = ffreqs[(nw1-nw2)/2];
            double v_R_s_low = ffreqs[(nw1-nw3)/2];

            vector<double> upperR_v1{v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_R_t_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v_K2_s_up,v_R_s_up};
            sort(upperR_v1.begin(),upperR_v1.end());
            vector<double> lowerR_v1{v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_R_t_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v_K2_s_low,v_R_s_low};
            sort(lowerR_v1.begin(),lowerR_v1.end());

            vector<double> upperR_v2{v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_R_t_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_s_up,v_R_s_up};
            sort(upperR_v2.begin(),upperR_v2.end());
            vector<double> lowerR_v2{v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_R_t_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_s_low,v_R_s_low};
            sort(lowerR_v2.begin(),lowerR_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,s,wlimit,w2,'s',1,h);
            vert2_const = vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,wlimit,'s',2,h);

            for(int i=upperR_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,upperR_v2[i],'s',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,s,upperR_v2[i],w2,'s',1,h))>0 ){limit_up2 = upperR_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperR_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,upperR_v1[i],w2,'s',1,h)-vert1_const)* vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperR_v1[i],'s',2,h))>0 ){limit_up1 = upperR_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            for(int i=0; i< lowerR_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerR_v1[i],w2,'s',1,h)-vert1_const)*vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerR_v1[i],'s',2,h))>0){limit_low1 = lowerR_v1[i];break;};
                if(i==lowerR_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerR_v2.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerR_v2[i],'s',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerR_v2[i],w2,'s',1,h))>0){limit_low2 = lowerR_v2[i];break;};
                if(i==lowerR_v2.size()-1){mode2_low=0;};
            };


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};}

        else if(h=='L'){//in this case, only Gamma (vert2) sets the limits

            double v_K2_s_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_s_up = ffreqs[(nw1+nw3)/2-1];
            double v_K2_s_low = ffreqs[(nw1-nw2)/2];
            double v_R_s_low = ffreqs[(nw1-nw3)/2];

            //conditions from vertex 2 (unprimed):
            double v2_K1_t_up = bfreqs[nw1-1]+w1;
            double v2_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw+nw2)/2-1]-w1-s},compare);
            double v2_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw+nw2)/2-1]-w1+s},compare);
            double v2_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw+nw3)/2-1]-w1+s,2*ffreqs[(nw+nw3)/2-1]-w1-s},compare);
            double v2_K1_u_up = bfreqs[nw1-1]-w1;
            double v2_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1+s},compare);
            double v2_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw+nw2)/2-1]+w1-s},compare);
            double v2_R_u_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw+nw3)/2-1]+w1+s,2*ffreqs[(nw+nw3)/2-1]+w1-s},compare);


            //conditions from vertex 2 (unprimed):
            double v2_K1_t_low =  bfreqs[0]+w1;
            double v2_K2_t_low = max({bfreqs[(nw-+nw2)/2]+w1,2*ffreqs[(nw-nw2)/2]-w1-s},compare);
            double v2_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw-nw2)/2]-w1+s},compare);
            double v2_R_t_low = max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw-nw3)/2]-w1+s,2*ffreqs[(nw-nw3)/2]-w1-s},compare);
            double v2_K1_u_low = bfreqs[0]-w1;
            double v2_K2_u_low =  max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1+s},compare);
            double v2_K2b_u_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw-nw2)/2]+w1-s},compare);
            double v2_R_u_low = max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw-nw3)/2]+w1+s,2*ffreqs[(nw-nw3)/2]+w1-s},compare);

            vector<double> upperK2_v1{v_K2_s_up};
            sort(upperK2_v1.begin(),upperK2_v1.end());
            vector<double> lowerK2_v1{v_K2_s_low};
            sort(lowerK2_v1.begin(),lowerK2_v1.end());

            vector<double> upperK2_v2{v2_K1_t_up,v2_K2_t_up,v2_K2b_t_up,v2_R_t_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_s_up,v_R_s_up};
            sort(upperK2_v2.begin(),upperK2_v2.end());
            vector<double> lowerK2_v2{v2_K1_t_low,v2_K2_t_low,v2_K2b_t_low,v2_R_t_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_s_low,v_R_s_low};
            sort(lowerK2_v2.begin(),lowerK2_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,s,wlimit,wlimit,'s',1,h);
            vert2_const = vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,wlimit,'s',2,h);


            for(int i=upperK2_v2.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK2_v2[i],'s',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK2_v2[i],w2,'s',1,h))>0 ){limit_up2 = upperK2_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK2_v1.size()-1; i>-1 ; i--){

                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK2_v1[i],w2,'s',1,h)-vert1_const)* vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK2_v1[i],'s',2,h))>0 ){limit_up1 = upperK2_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};

            for(int i=0; i< lowerK2_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK2_v1[i],w2,'s',1,h)-vert1_const)*vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK2_v1[i],'s',2,h))>0){limit_low1 = lowerK2_v1[i];break;};
                if(i==lowerK2_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK2_v2.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK2_v2[i],'s',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK2_v2[i],w2,'s',1,h))>0){limit_low2 = lowerK2_v2[i];break;};
                if(i==lowerK2_v2.size()-1){mode2_low=0;};
            };

            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};


            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};

        }

        else if(h=='M'){//in this case, only Gamma' (vert1) sets the limits

            //conditions from vertex 1 (primed):
            double v1_K1_t_up = bfreqs[nw1-1]+w2;
            double v1_K2_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-s},compare);//choose the minimum since beyond this value, this vertex cannot contribute any more
            double v1_K2b_t_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+s},compare);
            double v1_R_t_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2-s,2*ffreqs[(nw+nw3)/2-1]-w2+s},compare);
            double v1_K1_u_up= bfreqs[nw1-1]-w2;
            double v1_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-s},compare);
            double v1_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+s},compare);
            double v1_R_u_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2-s,2*ffreqs[(nw+nw3)/2-1]+w2+s},compare);

            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_t_low = bfreqs[0]+w2;
            double v1_K2_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-s},compare);//choose the minimum since bezond this value, this vertex cannot contribute any more
            double v1_K2b_t_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+s},compare);
            double v1_R_t_low = max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2-s,2*ffreqs[(nw-nw3)/2]-w2+s},compare);
            double v1_K1_u_low = bfreqs[0]-w2;
            double v1_K2_u_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-s},compare);
            double v1_K2b_u_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+s},compare);
            double v1_R_u_low = max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2-s,2*ffreqs[(nw-nw3)/2]+w2+s},compare);

            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_s_low = ffreqs[(nw1-nw2)/2];
            double v_R_s_low = ffreqs[(nw1-nw3)/2];
            //conditions that is equal for both vertices (s-channel contributions). Note that K_{s,1} is constant
            double v_K2_s_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_s_up = ffreqs[(nw1+nw3)/2-1];

            vector<double> upperK2b_v1{v1_K1_t_up,v1_K2_t_up,v1_K2b_t_up,v1_R_t_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v_K2_s_up,v_R_s_up};
            sort(upperK2b_v1.begin(),upperK2b_v1.end());
            vector<double> lowerK2b_v1{v1_K1_t_low,v1_K2_t_low,v1_K2b_t_low,v1_R_t_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v_K2_s_low,v_R_s_low};
            sort(lowerK2b_v1.begin(),lowerK2b_v1.end());

            vector<double> upperK2b_v2{v_K2_s_up};
            sort(upperK2b_v2.begin(),upperK2b_v2.end());
            vector<double> lowerK2b_v2{v_K2_s_low};
            sort(lowerK2b_v2.begin(),lowerK2b_v2.end());


            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,s,wlimit,w2,'s',1,h);
            vert2_const =vert1.vvalsmooth(red_side,map2,d,e,f,s,wlimit,wlimit,'s',2,h);

            for(int i=upperK2b_v2.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK2b_v2[i],'s',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK2b_v2[i],w2,'s',1,h))>0 ){limit_up2 = upperK2b_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK2b_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK2b_v1[i],w2,'s',1,h)-vert1_const)* vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK2b_v1[i],'s',2,h))>0 ){limit_up1 = upperK2b_v1[i];break;};
                if(i==0){mode1_up=0;};
            };


            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            for(int i=0; i< lowerK2b_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK2b_v1[i],w2,'s',1,h)-vert1_const)*vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK2b_v1[i],'s',2,h))>0){limit_low1 = lowerK2b_v1[i];break;};
                if(i==lowerK2b_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK2b_v2.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK2b_v2[i],'s',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK2b_v2[i],w2,'s',1,h))>0){limit_low2 = lowerK2b_v2[i];break;};
                if(i==lowerK2b_v2.size()-1){mode2_low=0;};
            };

            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};

        }

        else if(h=='K'){ //in this case the only contribution to the non-constat part comes from K2,t

            double v_K2_s_low = ffreqs[(nw1-nw2)/2];

            double v_K2_s_up = ffreqs[(nw1+nw2)/2-1];


            vector<double> upperK1_v1{v_K2_s_up};
            sort(upperK1_v1.begin(),upperK1_v1.end());
            vector<double> lowerK1_v1{v_K2_s_low};
            sort(lowerK1_v1.begin(),lowerK1_v1.end());

            vector<double> upperK1_v2{v_K2_s_up};
            sort(upperK1_v2.begin(),upperK1_v2.end());
            vector<double> lowerK1_v2{v_K2_s_low};
            sort(lowerK1_v2.begin(),lowerK1_v2.end());


            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,s,wlimit,wlimit,'s',1,h);
            vert2_const = vert1.vvalsmooth(red_side,map2,d,e,f,s,wlimit,wlimit,'s',2,h);



            for(int i=upperK1_v2.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK1_v2[i],'s',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK1_v2[i],w2,'s',1,h))>0 ){limit_up2 = upperK1_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK1_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,upperK1_v1[i],w2,'s',1,h)-vert1_const)* vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,upperK1_v1[i],'s',2,h))>0 ){limit_up1 = upperK1_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            for(int i=0; i< lowerK1_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK1_v1[i],w2,'s',1,h)-vert1_const)*vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK1_v1[i],'s',2,h))>0){limit_low1 = lowerK1_v1[i];break;};
                if(i==lowerK1_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK1_v2.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map2,d,e,f,s,w1,lowerK1_v2[i],'s',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,s,lowerK1_v2[i],w2,'s',1,h))>0){limit_low2 = lowerK1_v2[i];break;};
                if(i==lowerK1_v2.size()-1){mode2_low=0;};
            };
            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};


        };

        double result_re=0, error_re=0;
        double resultgfl_re=0, errorgfl_re=0;
        double resultgfu_re=0, errorgfu_re=0;

        struct pbubble_params<T1,T2> params = {red_side,map1,map2,Lambda,vert1, vert2,p1,p2,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};
        gsl_function F;

        if(mode==0 && abs(vert1_const * vert2_const)>0){

            if(REG==1 && p1=='g' && p2 =='g'){
                F.function = &pgreensfunc_re<T1,T2>;
                F.params = &params;

                //even integrand --> sufficient to intergrate one side and multiply result by two
                double upper = abs(s/2) + Lambda;


                F.function = &pgreensfunc_re<T1,T2>;
                gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                      w, &resultgfu_re, &errorgfu_re);
                double compl_real = 2. * resultgfu_re *vert1_const * vert2_const;

                B = compl_real;

            }


            else if(REG==1 && ((p1=='k' &&p2=='g') ||(p1=='g'&& p2 =='k'))){
                double lower = -abs(s/2)-Lambda;
                double upper = abs(s/2) + Lambda;
                char p1_new,p2_new;
                double dse_cutoff_up=0, dse_cutoff_low=0;//cutoff due to finite number of saved value in differentiated self energy
                double singsc=0;




                if(p1=='k'){
                    singsc =  real(1./2 * 1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,-Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,-Lambda-s/2,'s',2,h)  *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,s+Lambda,se,dse,p2)) ;
                    singsc += real(1./2 *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,s,Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,Lambda-s/2,'s',2,h)   *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,s-Lambda, se,dse, p2));
                    p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]-s/2;dse_cutoff_low= ffreqs[0]-s/2;}

                else if(p2=='k'){
                    singsc =   real(1./2 * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2+Lambda,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2+Lambda,'s',2,h)  *  propag(Lambda,s+Lambda,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
                    singsc += real(1./2 *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2-Lambda,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2-Lambda,'s',2,h)  *  propag(Lambda,s-Lambda,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
                    p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]+s/2;dse_cutoff_low= ffreqs[0]+s/2;};

                struct pbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(dse_cutoff_low <lower && dse_cutoff_up <=upper){//if highest saved dse-value if lower than upper heavyside-cutoff

                    if(dse_cutoff_up>lower){dse_cutoff_up=lower;};

                    F.function = &pgreensfunc_re<T1,T2>;
                    F.params = &params_mod;

                    F.function = &pgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    double compl_real = resultgfl_re *vert1_const * vert2_const;
                    B= singsc+ compl_real ;

                }

                else if(dse_cutoff_low >=lower && dse_cutoff_up > upper){//if lowest saved dse-value if higher than lower heavyside-cutoff

                    if(dse_cutoff_low<upper){dse_cutoff_low=upper;};

                    F.function = &pgreensfunc_re<T1,T2>;
                    F.params = &params_mod;

                    F.function = &pgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                    double compl_real = resultgfl_re *vert1_const * vert2_const;
                    B= singsc + compl_real;
                }
                else if (dse_cutoff_low<lower && dse_cutoff_up > upper){
                    F.function = &pgreensfunc_re<T1,T2>;
                    F.params = &params_mod;


                    F.function = &pgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,lower,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                    F.function = &pgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,upper,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfu_re, &errorgfu_re);
                    double compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    B= singsc+ compl_real ;
                };}


            else if(REG==2 && (p1=='g' &&p2=='g')){


                F.params = &params;
                F.function = &pgreensfunc_re<T1,T2>;
                gsl_integration_qagiu(&F,0,abs_error_bare,rel_error_bare,1500,
                                      w, &result_re, &errorgfu_re);

                double compl_real = 2*result_re *vert1_const * vert2_const;
                B=compl_real;

            }
            else if(REG==2 && ((p1=='g' &&p2=='s') ||(p1=='s'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially
                double bound_low=0, bound_up=0;

                if(p1=='s'){bound_low = -7*Lambda-s/2; bound_up = 7*Lambda-s/2;}
                else if(p2=='s'){bound_low = -7*Lambda+s/2; bound_up = 7*Lambda+s/2;}
                F.params = &params;
                F.function = &pgreensfunc_re<T1,T2>;
                gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                    w, &result_re, &error_re);
                double compl_real = result_re *vert1_const * vert2_const;
                B=compl_real;
            }

            else if(REG==2 && ((p1=='g' &&p2=='k') ||(p1=='k'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially
                double bound_low=0, bound_up=0, dse_cutoff_low=0, dse_cutoff_up=0;
                char p1_new, p2_new, p1_singsc,p2_singsc;
                if(p1=='k'){
                    bound_low = -7*Lambda-s/2; bound_up = 7*Lambda-s/2;
                    dse_cutoff_low = ffreqs[0]-s/2; dse_cutoff_up = ffreqs[nw-1]-s/2;
                    p1_new = 'e';
                    p1_singsc = 's';
                    p2_new = 'g';
                    p2_singsc = 'g';
                }
                else if(p2=='k'){
                    bound_low = -7*Lambda+s/2; bound_up = 7*Lambda+s/2;
                    dse_cutoff_low = ffreqs[0]+s/2; dse_cutoff_up = ffreqs[nw-1] +s/2;
                    p1_new = 'g';
                    p1_singsc = 'g';
                    p2_new = 'e';
                    p2_singsc = 's';
                }

                //first compute single scale contribution
                struct pbubble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_singsc,p2_singsc,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};
                F.params = &params_singsc;
                F.function = &pgreensfunc_re<T1,T2>;//integration of single scale
                gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);

                //add katanin extension
                struct pbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                F.params = &params_mod;//add integration of katanin extension
                F.function = &pgreensfunc_re<T1,T2>;
                double result_re2=0,error_re2=0;
                gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error, rel_error,1500,2,
                                    w, &result_re2, &error_re2);

                double compl_real = (result_re + result_re2 ) *vert1_const * vert2_const;

                B=compl_real;

            };
        }

        else if(mode==1){


            if(REG ==1 && p1=='g' && p2=='g'){
                double lhs=0, rhs=0;//result on left and right side of sharp cutoff
                double lower = -abs(s/2)-Lambda;
                double upper = abs(s/2) + Lambda;

                if(lower <= limit_low){//left side of sharp cutoff

                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,lower,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);

                        double compl_real = resultgfl_re  *vert1_const * vert2_const;

                        lhs=compl_real;

                    };}
                else if(lower <= limit_up){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);
                        compl_real = resultgfl_re *vert1_const * vert2_const;

                    };

                    F.params = &params;
                    F.function = &pbubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,lower,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    lhs=result_re + compl_real;

                }

                else if(lower > limit_up){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);

                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,lower,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                    };
                    F.params = &params;
                    F.function = &pbubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    lhs=result_re + compl_real;

                };

                if(upper >= limit_up){//right side of sharp cutoff

                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                        double compl_real = resultgfu_re *vert1_const * vert2_const;
                        rhs=compl_real;
                    };}
                else if(upper >= limit_low){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                        compl_real =resultgfu_re  *vert1_const * vert2_const;
                    };
                    F.params = &params;
                    F.function = &pbubble_re<T1,T2>;

                    gsl_integration_qag(&F,upper,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    rhs=result_re + compl_real;

                }

                else if (upper < limit_low){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,upper,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                    };
                    F.params = &params;
                    F.function = &pbubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    rhs = result_re + compl_real;
                };

                B =rhs + lhs;
            }

            else if(REG==1 && ((p1=='k' &&p2=='g') ||(p1=='g'&& p2 =='k'))){

                double lower = -abs(s/2)-Lambda;
                double upper = abs(s/2) + Lambda;
                char p1_new,p2_new;
                double dse_cutoff_up=0, dse_cutoff_low=0;//cutoff due to finite number of saved value in differentiated self energy
                double singsc=0;


                if(p1=='k'){
                    singsc =  real(1./2 * 1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,-Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,-Lambda-s/2,'s',2,h)  *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,s+Lambda,se,dse,p2)) ;
                    singsc += real(1./2 *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,s,Lambda-s/2,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,Lambda-s/2,'s',2,h)   *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,s-Lambda, se,dse, p2));
                    p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]-s/2;dse_cutoff_low= ffreqs[0]-s/2;}

                else if(p2=='k'){
                    singsc =   real(1./2 * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2+Lambda,w2,'s',1,h) *  vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2+Lambda,'s',2,h)  *  propag(Lambda,s+Lambda,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
                    singsc += real(1./2 *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,s,s/2-Lambda,w2,'s',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,s,w1,s/2-Lambda,'s',2,h)  *  propag(Lambda,s-Lambda,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
                    p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]+s/2;dse_cutoff_low= ffreqs[0]+s/2;};


                double lhs=0, rhs=0;//result on left and right side of sharp cutoff

                struct pbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(lower <= limit_low && abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){//left side of sharp cutoff

                    double up_eff=0;
                    if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                    F.params = &params_mod;
                    F.function = &pgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low ,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);
                    double compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    lhs=compl_real;

                }
                else if(lower <= limit_up ){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };
                    if(dse_cutoff_low <lower){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                        F.params = &params_mod;
                        F.function = &pbubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                        lhs=result_re + compl_real;
                    }

                    else if(lower > limit_up){

                        double compl_real=0;
                        if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                            double up_eff=0;
                            if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                            F.params = &params_mod;
                            F.function = &pgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfl_re, &errorgfl_re);
                            compl_real = resultgfl_re  *vert1_const * vert2_const;
                        };


                        if(dse_cutoff_low <limit_up){
                            double low_eff=0;
                            if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                            double up_eff=0;
                            if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                            F.params = &params_mod;
                            F.function = &pbubble_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                                w, &result_re, &error_re);
                        };


                        if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){
                            double low_eff=0;
                            if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                            double up_eff=0;
                            if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                            F.params = &params_mod;
                            F.function = &pgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfl_re, &errorgfl_re);
                            compl_real += resultgfl_re  *vert1_const * vert2_const;

                        };

                        lhs=result_re + compl_real;
                    };

                    if(upper >= limit_up && abs(vert1_const * vert2_const)>0 && dse_cutoff_up >upper){//right side of sharp cutoff

                        double low_eff=0;
                        if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                        F.params = &params_mod;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        double compl_real = (resultgfu_re  )*vert1_const * vert2_const;

                        rhs=compl_real;

                    }
                    else if(upper >= limit_low ){
                        double compl_real=0;
                        if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up >limit_up){
                            double low_eff=0;
                            if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                            F.params = &params_mod;
                            F.function = &pgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfu_re, &errorgfu_re);
                            compl_real =(resultgfu_re  )*vert1_const * vert2_const;
                        };
                        if(dse_cutoff_up > upper){
                            double low_eff=0;
                            if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                            double up_eff=0;
                            if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                            F.params = &params_mod;
                            F.function = &pbubble_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                                w, &result_re, &error_re);

                        };

                        rhs=result_re + compl_real;

                    }

                    else if(upper < limit_low){

                        double compl_real=0;
                        if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > upper){
                            double low_eff=0;
                            if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                            double up_eff=0;
                            if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                            F.params = &params_mod;

                            F.function = &pgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfu_re, &errorgfu_re);
                            compl_real = (resultgfu_re  )*vert1_const * vert2_const;
                        };


                        if(dse_cutoff_up >limit_low){
                            double low_eff=0;
                            if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                            double up_eff=0;
                            if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                            F.params = &params_mod;
                            F.function = &pbubble_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                                w, &result_re, &error_re);
                        };


                        if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > limit_up){
                            double low_eff=0;
                            if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};

                            F.params = &params_mod;
                            F.function = &pgreensfunc_re<T1,T2>;
                            gsl_integration_qag(&F,low_eff,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                                w, &resultgfu_re, &errorgfu_re);
                            compl_real += (resultgfu_re  )*vert1_const * vert2_const;

                        };

                        rhs=result_re + compl_real;

                    };


                    B = rhs + lhs + singsc;
                };}


            else if(REG==2 && p1=='g' && p2=='g'){

                F.function = &pbubble_re<T1,T2>;
                F.params = &params;

                gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);
                //integration of greens function outside of w-dependent interval

                double compl_real=0;

                if(abs(vert1_const * vert2_const)>0){
                    F.function = &pgreensfunc_re<T1,T2>;
                    gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfl_re, &errorgfl_re);

                    F.function = &pgreensfunc_re<T1,T2>;
                    gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfu_re, &errorgfu_re);

                    compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                };
                B=result_re + compl_real ;
            }

            else if(REG==2 && ((p1=='g' && p2=='s') || (p1=='s' && p2=='g'))){

                double bound_low=0, bound_up=0;

                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;

                if(p1=='s'){bound_low = -7*Lambda-s/2; bound_up = 7*Lambda-s/2;}
                else if(p2=='s'){bound_low = -7*Lambda+s/2; bound_up = 7*Lambda+s/2;}

                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                    };




                    B=result_re + compl_real ;
                }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real = (resultgfl_re )*vert1_const * vert2_const;
                    };


                    B = result_re + compl_real ;
                }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){

                        F.params = &params;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfu_re )*vert1_const * vert2_const;
                    };



                    B=result_re + compl_real ;
                }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    B=compl_real ;
                }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    B=result_re  ;
                }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real =(resultgfu_re )*vert1_const * vert2_const;
                    };

                    B=compl_real ;
                }

            }

            else if(REG==2 && ((p1=='g' && p2=='k') || (p1=='k' && p2=='g'))){

                double singsc=0;
                double bound_low=0, bound_up=0;
                //FIRST COMPUTE THE SINGLE SCALE CONTRIBUTION
                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;
                char p1_new,p2_new;
                if(p1=='k'){
                    bound_low = -7*Lambda-s/2; bound_up = 7*Lambda-s/2;
                    p1_new = 's'; p2_new = 'g';}
                else if(p2=='k'){bound_low = -7*Lambda+s/2; bound_up = 7*Lambda+s/2;
                    p1_new = 'g'; p2_new = 's';
                }


                struct pbubble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                    };




                    singsc=result_re + compl_real;
                }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re)*vert1_const * vert2_const;

                    };


                    singsc=result_re + compl_real ;
                }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){

                        F.params = &params_singsc;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfu_re )*vert1_const * vert2_const;
                    };



                    singsc=result_re + compl_real ;
                }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real = (resultgfl_re )*vert1_const * vert2_const;
                    };


                    singsc= compl_real ;
                }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    singsc=result_re  ;
                }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfu_re)*vert1_const * vert2_const;

                    };

                    singsc= compl_real ;
                }
                //ADD THE CONTRIBUTION FROM THE KATANON EXTENSION:
                double kat=0;

                //reset m-values
                m1=0;//turns on and off the integration in the three parts
                m2=0;
                m3=0;


                if(p1=='k'){
                    bound_low = ffreqs[0]-s/2; bound_up = ffreqs[nw-1]-s/2;//the bounds are the boundaries of the dse-grid
                    p1_new = 'e';
                    p2_new = 'g';}
                else if(p2=='k'){
                    bound_low = ffreqs[0]+s/2; bound_up = ffreqs[nw-1]+s/2;
                    p1_new = 'g';
                    p2_new = 'e';}

                struct pbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  s, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };




                    kat=result_re + compl_real ;
                }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real = (resultgfl_re )*vert1_const * vert2_const;

                    };


                    kat=(result_re + compl_real );
                }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){


                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    };



                    kat=(result_re + compl_real );
                }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    //    double compl_imag;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real =(resultgfl_re )*vert1_const * vert2_const;

                    };


                    kat=( compl_real );
                }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &pbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    kat=result_re  ;
                }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &pgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    };

                    kat= compl_real ;
                }
                B=kat + singsc;

            };


        };

    };
    if(abs(B)<1e-20){B=0;};
    return B;
}



//t-bubble:
template<typename  T1,typename  T2>
struct tbubble_params{
    int red_side;
    int map1;
    int map2;
    double Lambda;
    T1& vert1;
    T2& vert2;
    char ptype1;
    char ptype2;
    SelfEnergy<comp>& selfenergy;
    SelfEnergy<comp>& diffselfenergy;

    int a; int b; int c;

    int d; int e; int f;

    double t;double w1; double w2;
    char h;

};

template<typename T1,typename T2>
double tbubble_re(double w, void * p) {
    struct tbubble_params<T1,T2> * params
            = static_cast< struct tbubble_params<T1,T2> *>(p);
    int red_side = (params->red_side);
    int map1 = (params->map1);
    int map2 = (params->map2);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
    SelfEnergy<comp>& selfenergy = (params->selfenergy);
    SelfEnergy<comp>& diffselfenergy = (params->diffselfenergy);
    double Lambda = (params->Lambda);
    T1& vert1 = (params-> vert1);
    T2& vert2 = (params -> vert2);
    int a  = (params-> a);
    int b = (params->b);
    int c = (params->c);

    int d =(params->d);
    int  e= (params->e);
    int f = (params->f);

    double t = (params->t);
    double w1 = (params->w1);
    double w2 = (params->w2);
    char h = (params -> h);

    double val = real(-1.*  vert1.vvalsmooth(red_side,map1,a,b,c,t,w,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,w,'t',2,h) *  propag(Lambda,w-t/2,selfenergy,diffselfenergy,ptype1) * propag(Lambda,w+t/2,selfenergy,diffselfenergy,ptype2) );

    return (1./(2*pi)*val);
}

template<typename T1,typename T2>
double tgreensfunc_re(double w, void * p){//returns the value of all constituents of a bubble function at a specific value of the integration parameter WITHOUT contribution of bare vertices!
    struct tbubble_params<T1,T2> * params
            = static_cast< struct tbubble_params<T1,T2> *>(p);
    char ptype1 = (params ->ptype1);
    char ptype2 = (params ->ptype2);
    SelfEnergy<comp>& selfenergy = (params->selfenergy);
    SelfEnergy<comp>& diffselfenergy = (params->diffselfenergy);
    double Lambda = (params->Lambda);

    double t = (params->t);


    double val = real( propag(Lambda,w-t/2,selfenergy,diffselfenergy,ptype1) * propag(Lambda,w+t/2,selfenergy,diffselfenergy,ptype2) );


    return (-1./(2.*pi)*val);//minus sign is from definition of t-bubble
}

template<typename T1,typename T2>
double tbubble(int red_side, int map1, int map2,gsl_integration_workspace* w,double Lambda, T1& vert1,int a, int b, int c, T2& vert2,int d, int e, int f,char p1, char p2, SelfEnergy<comp>& se, SelfEnergy<comp>& dse, double t,double w1, double w2, char h){//bubble in s-channel: specified by s: bos. frequency and two external ferm freqs.

    double abs_error=1e-2, rel_error=1e-4;

    double abs_error_bare=1e-2,rel_error_bare=1e-4;

    double B=0;


    if((REG ==1 && p1 =='s' && p2 =='g' ) ){//if first propagator is single scale propagator, proportional to delta peak: no integration
        B = real( -1. *1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda+t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda+t/2,'t',2,h) * propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+t,se,dse,p2)) ;
        B += real(-1. *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda+t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda+t/2,'t',2,h) * propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,Lambda+t,se,dse,p2)) ;
    }
    else if((REG ==1 && p1=='g' && p2 =='s') ){//if second propagator is single scale propagator, proportional to delta peak: no integration

        B = real(-1. * 1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda-t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda-t/2,'t',2,h) * propag(Lambda,-Lambda-t,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
        B += real(-1. *1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda-t/2,w2,'t',1,h)  *   vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda-t/2,'t',2,h) * propag(Lambda,Lambda-t,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
    }

    else{

        double limit_low=0,limit_low1=0,limit_low2=0;
        double limit_up=0,limit_up1=0,limit_up2=0;

        double vert1_const,vert2_const;





        int mode =1;

        if(h == 'R'){

            //upper bound:

            //conditions from vertex 1 (primed):
            double v1_K1_u_up = bfreqs[nw1-1]+w2;
            double v1_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+t},compare);
            double v1_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-t},compare);
            double v1_R_u_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2+t,2*ffreqs[(nw+nw3)/2-1]-w2-t},compare);
            double v1_K1_s_up= bfreqs[nw1-1]-w2;
            double v1_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+t},compare);
            double v1_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-t},compare);
            double v1_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2+t,2*ffreqs[(nw+nw3)/2-1]+w2-t},compare);

            //conditions from vertex 2 (unprimed):
            double v2_K1_u_up = bfreqs[nw1-1]+w1;
            double v2_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1+t},compare);
            double v2_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1-t},compare);
            double v2_R_u_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw1+nw3)/2-1]-w1+t,2*ffreqs[(nw1+nw3)/2-1]-w1-t},compare);
            double v2_K1_s_up = bfreqs[nw1-1]-w1;
            double v2_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw1+nw2)/2-1]+w1-t},compare);
            double v2_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw1+nw2)/2-1]+w1+t},compare);
            double v2_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw1+nw3)/2-1]+w1+t,2*ffreqs[(nw1+nw3)/2-1]+w1-t},compare);
            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_t_up = ffreqs[(nw1+nw3)/2-1];



            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_u_low = bfreqs[0]+w2;
            double v1_K2_u_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+t},compare);
            double v1_K2b_u_low =  max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-t},compare);
            double v1_R_u_low =  max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2+t,2*ffreqs[(nw-nw3)/2]-w2-t},compare);
            double v1_K1_s_low = bfreqs[0]-w2;
            double v1_K2_s_low =  max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+t},compare);
            double v1_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-t},compare);
            double v1_R_s_low =  max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2+t,2*ffreqs[(nw-nw3)/2]+w2-t},compare);

            //conditions from vertex 2 (unprimed):
            double v2_K1_u_low = bfreqs[0]+w1;
            double v2_K2_u_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1+t},compare);
            double v2_K2b_u_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1-t},compare);
            double v2_R_u_low = max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw1-nw3)/2]-w1+t,2*ffreqs[(nw1-nw3)/2]-w1-t},compare);
            double v2_K1_s_low =  bfreqs[0]-w1;
            double v2_K2_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw1-nw2)/2]+w1-t},compare);
            double v2_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw1-nw2)/2]+w1+t},compare);
            double v2_R_s_low =max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw1-nw3)/2]+w1+t,2*ffreqs[(nw1-nw3)/2]+w1-t},compare);
            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_low = ffreqs[(nw1-nw2)/2];
            double v_R_t_low = ffreqs[(nw1-nw3)/2];




            vector<double> upperR_v1{v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_R_s_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v_K2_t_up,v_R_t_up};
            sort(upperR_v1.begin(),upperR_v1.end());
            vector<double> lowerR_v1{v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_R_s_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v_K2_t_low,v_R_t_low};
            sort(lowerR_v1.begin(),lowerR_v1.end());

            vector<double> upperR_v2{v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_R_s_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_t_up,v_R_t_up};
            sort(upperR_v2.begin(),upperR_v2.end());
            vector<double> lowerR_v2{v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_R_s_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_t_low,v_R_t_low};
            sort(lowerR_v2.begin(),lowerR_v2.end());


            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,t,wlimit,w2,'t',1,h);
            vert2_const =vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,wlimit,'t',2,h);

            for(int i=upperR_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperR_v2[i],'t',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,t,upperR_v2[i],w2,'t',1,h))>0 ){limit_up2 = upperR_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperR_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,upperR_v1[i],w2,'t',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperR_v1[i],'t',2,h))>0 ){limit_up1 = upperR_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};



            for(int i=0; i< lowerR_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerR_v1[i],w2,'t',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerR_v1[i],'t',2,h))>0){limit_low1 = lowerR_v1[i];break;};
                if(i==lowerR_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerR_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerR_v2[i],'t',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerR_v2[i],w2,'t',1,h))>0){limit_low2 = lowerR_v2[i];break;};
                if(i==lowerR_v2.size()-1){mode2_low=0;};
            };


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};}

        else if(h == 'L'){ //constant vertices


            //conditions from vertex 2 (unprimed):
            double v2_K1_u_up = bfreqs[nw1-1]+w1;
            double v2_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1+t},compare);
            double v2_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]+w1,2*ffreqs[(nw1+nw2)/2-1]-w1-t},compare);
            double v2_R_u_up = min({bfreqs[(nw1+nw3)/2-1]+w1,2*ffreqs[(nw1+nw3)/2-1]-w1+t,2*ffreqs[(nw1+nw3)/2-1]-w1-t},compare);
            double v2_K1_s_up = bfreqs[nw1-1]-w1;
            double v2_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw1+nw2)/2-1]+w1-t},compare);
            double v2_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w1,2*ffreqs[(nw1+nw2)/2-1]+w1+t},compare);
            double v2_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w1,2*ffreqs[(nw1+nw3)/2-1]+w1+t,2*ffreqs[(nw1+nw3)/2-1]+w1-t},compare);

            //conditions from vertex 2 (unprimed):
            double v2_K1_u_low = bfreqs[0]+w1;
            double v2_K2_u_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1+t},compare);
            double v2_K2b_u_low = max({bfreqs[(nw1-nw2)/2]+w1,2*ffreqs[(nw1-nw2)/2]-w1-t},compare);
            double v2_R_u_low = max({bfreqs[(nw1-nw3)/2]+w1,2*ffreqs[(nw1-nw3)/2]-w1+t,2*ffreqs[(nw1-nw3)/2]-w1-t},compare);
            double v2_K1_s_low =  bfreqs[0]-w1;
            double v2_K2_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw1-nw2)/2]+w1-t},compare);
            double v2_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w1,2*ffreqs[(nw1-nw2)/2]+w1+t},compare);
            double v2_R_s_low =max({bfreqs[(nw1-nw3)/2]-w1,2*ffreqs[(nw1-nw3)/2]+w1+t,2*ffreqs[(nw1-nw3)/2]+w1-t},compare);



            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_t_up = ffreqs[(nw1+nw3)/2-1];
            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_low = ffreqs[(nw1-nw2)/2];
            double v_R_t_low = ffreqs[(nw1-nw3)/2];




            vector<double> upperK2_v1{v_K2_t_up};

            vector<double> lowerK2_v1{v_K2_t_low};


            vector<double> upperK2_v2{v2_K1_s_up,v2_K2_s_up,v2_K2b_s_up,v2_R_s_up,v2_K1_u_up,v2_K2_u_up,v2_K2b_u_up,v2_R_u_up,v_K2_t_up,v_R_t_up};
            sort(upperK2_v2.begin(),upperK2_v2.end());
            vector<double> lowerK2_v2{v2_K1_s_low,v2_K2_s_low,v2_K2b_s_low,v2_R_s_low,v2_K1_u_low,v2_K2_u_low,v2_K2b_u_low,v2_R_u_low,v_K2_t_low,v_R_t_low};
            sort(lowerK2_v2.begin(),lowerK2_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,t,wlimit,wlimit,'t',1,h);
            vert2_const = vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,wlimit,'t',2,h);

            for(int i=upperK2_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK2_v2[i],'t',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK2_v2[i],w2,'t',1,h))>0 ){limit_up2 = upperK2_v2[i];break;};
                if(i==0){mode2_up=0;};


            };



            if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK2_v1[0],w2,'t',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK2_v1[0],'t',2,h))>0 ){limit_up1 = upperK2_v1[0];}
            else{mode1_up=0;};

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK2_v1[0],w2,'t',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK2_v1[0],'t',2,h))>0){limit_low1 = lowerK2_v1[0];}
            else{mode1_low=0;};


            for(int i=0; i< lowerK2_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK2_v2[i],'t',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK2_v2[i],w2,'t',1,h))>0){limit_low2 = lowerK2_v2[i];break;};
                if(i==lowerK2_v2.size()-1){mode2_low=0;};
            };


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};


        }

        else if(h == 'M'){//in this case, only Gamma' (vert1) sets the limits

            //lower bound:
            //conditions from vertex 1 (primed):
            double v1_K1_u_low = bfreqs[0]+w2;
            double v1_K2_u_low = max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2+t},compare);
            double v1_K2b_u_low =  max({bfreqs[(nw1-nw2)/2]+w2,2*ffreqs[(nw-nw2)/2]-w2-t},compare);
            double v1_R_u_low =  max({bfreqs[(nw1-nw3)/2]+w2,2*ffreqs[(nw-nw3)/2]-w2+t,2*ffreqs[(nw-nw3)/2]-w2-t},compare);
            double v1_K1_s_low = bfreqs[0]-w2;
            double v1_K2_s_low =  max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2+t},compare);
            double v1_K2b_s_low = max({bfreqs[(nw1-nw2)/2]-w2,2*ffreqs[(nw-nw2)/2]+w2-t},compare);
            double v1_R_s_low =  max({bfreqs[(nw1-nw3)/2]-w2,2*ffreqs[(nw-nw3)/2]+w2+t,2*ffreqs[(nw-nw3)/2]+w2-t},compare);

            //conditions from vertex 1 (primed):
            double v1_K1_u_up = bfreqs[nw1-1]+w2;
            double v1_K2_u_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2+t},compare);
            double v1_K2b_u_up = min({bfreqs[(nw1+nw2)/2-1]+w2,2*ffreqs[(nw+nw2)/2-1]-w2-t},compare);
            double v1_R_u_up = min({bfreqs[(nw1+nw3)/2-1]+w2,2*ffreqs[(nw+nw3)/2-1]-w2+t,2*ffreqs[(nw+nw3)/2-1]-w2-t},compare);
            double v1_K1_s_up= bfreqs[nw1-1]-w2;
            double v1_K2_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2+t},compare);
            double v1_K2b_s_up = min({bfreqs[(nw1+nw2)/2-1]-w2,2*ffreqs[(nw+nw2)/2-1]+w2-t},compare);
            double v1_R_s_up = min({bfreqs[(nw1+nw3)/2-1]-w2,2*ffreqs[(nw+nw3)/2-1]+w2+t,2*ffreqs[(nw+nw3)/2-1]+w2-t},compare);


            double v_K2_t_up = ffreqs[(nw1+nw2)/2-1];
            double v_R_t_up = ffreqs[(nw1+nw3)/2-1];


            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_low = ffreqs[(nw1-nw2)/2];
            double v_R_t_low = ffreqs[(nw1-nw3)/2];

            vector<double> upperK2b_v1{v1_K1_s_up,v1_K2_s_up,v1_K2b_s_up,v1_R_s_up,v1_K1_u_up,v1_K2_u_up,v1_K2b_u_up,v1_R_u_up,v_K2_t_up,v_R_t_up};
            sort(upperK2b_v1.begin(),upperK2b_v1.end());
            vector<double> lowerK2b_v1{v1_K1_s_low,v1_K2_s_low,v1_K2b_s_low,v1_R_s_low,v1_K1_u_low,v1_K2_u_low,v1_K2b_u_low,v1_R_u_low,v_K2_t_low,v_R_t_low};
            sort(lowerK2b_v1.begin(),lowerK2b_v1.end());

            vector<double> upperK2b_v2{v_K2_t_up};
            sort(upperK2b_v2.begin(),upperK2b_v2.end());
            vector<double> lowerK2b_v2{v_K2_t_low};
            sort(lowerK2b_v2.begin(),lowerK2b_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;

            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,t,wlimit,w2,'t',1,h);
            vert2_const =vert2.vvalsmooth(red_side,map2,d,e,f,t,wlimit,wlimit,'t',2,h);

            for(int i=upperK2b_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK2b_v2[i],'t',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK2b_v2[i],w2,'t',1,h))>0 ){limit_up2 = upperK2b_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK2b_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK2b_v1[i],w2,'t',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK2b_v1[i],'t',2,h))>0 ){limit_up1 = upperK2b_v1[i];break;};
                if(i==0){mode1_up=0;};
            };


            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};

            for(int i=0; i< lowerK2b_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK2b_v1[i],w2,'t',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK2b_v1[i],'t',2,h))>0){limit_low1 = lowerK2b_v1[i];break;};
                if(i==lowerK2b_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK2b_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK2b_v2[i],'t',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK2b_v2[i],w2,'t',1,h))>0){limit_low2 = lowerK2b_v2[i];break;};
                if(i==lowerK2b_v2.size()-1){mode2_low=0;};
            };


            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};
            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};}

        else if(h == 'K'){ //in this case the only contribution to the non-constat part comes from K2,t

            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_low = ffreqs[(nw1-nw2)/2];

            //conditions that is equal for both vertices (t-channel contributions). Note that K_{t,1} is constant
            double v_K2_t_up = ffreqs[(nw1+nw2)/2-1];



            vector<double> upperK1_v1{v_K2_t_up};
            sort(upperK1_v1.begin(),upperK1_v1.end());
            vector<double> lowerK1_v1{v_K2_t_low};
            sort(lowerK1_v1.begin(),lowerK1_v1.end());

            vector<double> upperK1_v2{v_K2_t_up};
            sort(upperK1_v2.begin(),upperK1_v2.end());
            vector<double> lowerK1_v2{v_K2_t_low};
            sort(lowerK1_v2.begin(),lowerK1_v2.end());

            int mode1_up=1;int mode2_up=1;
            int mode1_low=1;int mode2_low=1;
            vert1_const = vert1.vvalsmooth(red_side,map1,a,b,c,t,wlimit,wlimit,'t',1,h);
            vert2_const = vert2.vvalsmooth(red_side,map2,d,e,f,t,wlimit,wlimit,'t',2,h);



            for(int i=upperK1_v2.size()-1; i>-1 ; i--){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK1_v2[i],'t',2,h)-vert2_const)* vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK1_v2[i],w2,'t',1,h))>0 ){limit_up2 = upperK1_v2[i];break;};
                if(i==0){mode2_up=0;};
            };

            for(int i=upperK1_v1.size()-1; i>-1 ; i--){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,upperK1_v1[i],w2,'t',1,h)-vert1_const)* vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,upperK1_v1[i],'t',2,h))>0 ){limit_up1 = upperK1_v1[i];break;};
                if(i==0){mode1_up=0;};
            };

            if(mode1_up == 1 && mode2_up ==1){
                if(limit_up1>limit_up2){limit_up=limit_up1;}
                else{limit_up=limit_up2;}}
            else if (mode1_up == 1 && mode2_up ==0){limit_up = limit_up1;}
            else if (mode1_up == 0 && mode2_up ==1){limit_up = limit_up2;};


            for(int i=0; i< lowerK1_v1.size() ; i++){
                if(abs((vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK1_v1[i],w2,'t',1,h)-vert1_const)*vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK1_v1[i],'t',2,h))>0){limit_low1 = lowerK1_v1[i];break;};
                if(i==lowerK1_v1.size()-1){mode1_low=0;};
            };

            for(int i=0; i< lowerK1_v2.size() ; i++){
                if(abs((vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,lowerK1_v2[i],'t',2,h)-vert2_const)*vert1.vvalsmooth(red_side,map1,a,b,c,t,lowerK1_v2[i],w2,'t',1,h))>0){limit_low2 = lowerK1_v2[i];break;};
                if(i==lowerK1_v2.size()-1){mode2_low=0;};
            };

            if(mode1_low == 1 && mode2_low ==1){
                if(limit_low1<limit_low2){limit_low=limit_low1;}
                else{limit_low=limit_low2;};}
            else if (mode1_low == 1 && mode2_low ==0){limit_low = limit_low1;}
            else if (mode1_low == 0 && mode2_low ==1){limit_low = limit_low2;};

            if(mode1_low ==0 && mode2_low ==0 && mode1_up ==0 && mode2_up ==0){mode =0;};



        };


        double result_re=0, error_re=0;
        double resultgfl_re=0, errorgfl_re=0;
        double resultgfu_re=0, errorgfu_re=0;

        struct tbubble_params<T1,T2> params = {red_side,map1,map2,Lambda,vert1, vert2,p1,p2,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};
        gsl_function F;





        if(mode==0 && abs(vert1_const * vert2_const)>0){

            if(REG==1 && p1=='g' && p2 =='g'){

                F.function = &tgreensfunc_re<T1,T2>;
                F.params = &params;


                double upper = abs(t/2) + Lambda;

//even integrand -> integrate only upper half and multiply by two


                F.function = &tgreensfunc_re<T1,T2>;
                gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                      w, &resultgfu_re, &errorgfu_re);

                double compl_real = 2*resultgfu_re *vert1_const * vert2_const;

                B= compl_real ;

            }
            else if(REG==1 && ((p1=='k' &&p2=='g') ||(p1=='g'&& p2 =='k'))){
                double lower = -abs(t/2)-Lambda;
                double upper = abs(t/2) + Lambda;



                char p1_new,p2_new;
                double dse_cutoff_up, dse_cutoff_low;//cutoff due to finite number of saved value in differentiated self energy
                double singsc;

                if(p1=='k'){
                    singsc = real( -1./(2*pi)* vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda+t/2,'t',2,h) *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+t,se,dse,p2)) ;
                    singsc += real(-1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda+t/2,'t',2,h) *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,Lambda+t,se,dse,p2)) ;
                    p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]+t/2;dse_cutoff_low= ffreqs[0]+t/2;}

                else if(p2=='k'){
                    singsc =  real(-1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda-t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda-t/2,'t',2,h) *  propag(Lambda,-Lambda-t,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
                    singsc += real(-1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda-t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda-t/2,'t',2,h) *  propag(Lambda,Lambda-t,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
                    p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]-t/2;dse_cutoff_low= ffreqs[0]-t/2;};

                struct tbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(dse_cutoff_low <lower && dse_cutoff_up <=upper){//if highest saved dse-value if lower than upper heavyside-cutoff

                    if(dse_cutoff_up>lower){dse_cutoff_up=lower;};


                    F.params = &params_mod;

                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    double compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    B=singsc+ compl_real ;
                }

                else if(dse_cutoff_low >=lower && dse_cutoff_up > upper){//if lowest saved dse-value if higher than lower heavyside-cutoff

                    if(dse_cutoff_low<upper){dse_cutoff_low=upper;};


                    F.params = &params_mod;

                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    double compl_real =(resultgfl_re  )*vert1_const * vert2_const;

                    B= singsc + compl_real ;
                }
                else if (dse_cutoff_low<lower && dse_cutoff_up > upper){

                    F.params = &params_mod;


                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low,lower,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,upper,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfu_re, &errorgfu_re);


                    double compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                    B= singsc ;//+ compl_real ;
                };}


            else if(REG==2 && (p1=='g' &&p2=='g')){
                F.params = &params;


                F.function = &tgreensfunc_re<T1,T2>;
                gsl_integration_qagiu(&F,0,abs_error_bare,rel_error_bare,1500,
                                      w, &result_re, &errorgfu_re);


                double compl_real = 2*result_re *vert1_const * vert2_const;

                B=compl_real;

            }
            else if(REG==2 && ((p1=='g' &&p2=='s') ||(p1=='s'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially
                double bound_low, bound_up;
                if(p1=='s'){bound_low = -7*Lambda+t/2; bound_up = 7*Lambda+t/2;}
                else if(p2=='s'){bound_low = -7*Lambda-t/2; bound_up = 7*Lambda-t/2;}
                F.params = &params;
                F.function = &tgreensfunc_re<T1,T2>;
                gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                    w, &result_re, &error_re);
                double compl_real = (result_re) *vert1_const * vert2_const;

                B=compl_real;

            }

            else if(REG==2 && ((p1=='g' &&p2=='k') ||(p1=='k'&& p2 =='g'))){//in this case, the integration interval must be very small around the peak for accurate numerical results since the single scale propag decays exponentially
                double bound_low, bound_up, dse_cutoff_low, dse_cutoff_up;
                char p1_new, p2_new, p1_singsc,p2_singsc;
                if(p1=='k'){
                    bound_low = -7*Lambda+t/2; bound_up = 7*Lambda+t/2;
                    dse_cutoff_low = ffreqs[0]+t/2; dse_cutoff_up = ffreqs[nw-1]+t/2;
                    p1_new = 'e';
                    p1_singsc = 's';
                    p2_new = 'g';
                    p2_singsc = 'g';
                }
                else if(p2=='k'){
                    bound_low = -7*Lambda-t/2; bound_up = 7*Lambda-t/2;
                    dse_cutoff_low = ffreqs[0]-t/2; dse_cutoff_up = ffreqs[nw-1] -t/2;
                    p1_new = 'g';
                    p1_singsc = 'g';
                    p2_new = 'e';
                    p2_singsc = 's';
                }

                //first compute single scale contribution
                struct tbubble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_singsc,p2_singsc,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};
                F.params = &params_singsc;
                F.function = &tgreensfunc_re<T1,T2>;//integration of single scale
                gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                    w, &result_re, &error_re);

                //add katanin extension
                struct tbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                F.params = &params_mod;//add integration of katanin extension
                F.function = &tgreensfunc_re<T1,T2>;
                double result_re2,error_re2,result_im2,error_im2;
                gsl_integration_qag(&F,dse_cutoff_low,dse_cutoff_up,abs_error_bare, rel_error_bare,1500,2,
                                    w, &result_re2, &error_re2);

                double compl_real =(result_re + result_re2 ) *vert1_const * vert2_const;

                B=compl_real;

            };
        }
        else if(mode==1){


            if(REG==1 && p1=='g' && p2=='g'){
                double lhs=0, rhs=0;//result on left and right side of sharp cutoff
                double lower = -abs(t/2)-Lambda;
                double upper = abs(t/2) + Lambda;

                if(lower <= limit_low){//left side of sharp cutoff

                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,lower,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);
                        double compl_real = (resultgfl_re )*vert1_const * vert2_const;
                        lhs=compl_real;
                    };}
                else if(lower <= limit_up){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);
                        compl_real = (resultgfl_re  )*vert1_const * vert2_const;
                    };

                    F.params = &params;
                    F.function = &tbubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,lower,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    lhs= result_re + compl_real;
                }

                else if(lower > limit_up){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfl_re, &errorgfl_re);

                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,lower,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };
                    F.params = &params;
                    F.function = &tbubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    lhs= result_re + compl_real;


                };

                if(upper >= limit_up){//right side of sharp cutoff

                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,upper,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                        double compl_real =(resultgfu_re  )*vert1_const * vert2_const;

                        rhs=(compl_real);

                    };}
                else if(upper >= limit_low){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfu_re  )*vert1_const * vert2_const;
                    };
                    F.params = &params;
                    F.function = &tbubble_re<T1,T2>;

                    gsl_integration_qag(&F,upper,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    rhs=result_re + compl_real;

                }

                else if (upper < limit_low){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,upper,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                              w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                    };
                    F.params = &params;
                    F.function = &tbubble_re<T1,T2>;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    rhs=result_re + compl_real;

                };
                B= rhs +lhs;
            }
            else if((REG==1 && p1=='k' && p2=='g') || (REG==1 && p1=='g' && p2=='k') ){


                double lower = -abs(t/2)-Lambda;
                double upper = abs(t/2) + Lambda;


                char p1_new,p2_new;
                double dse_cutoff_up=0, dse_cutoff_low=0;//cutoff due to finite number of saved value in differentiated self energy
                double singsc=0;

                if(p1=='k'){

                    singsc = real(- 1./(2*pi)*  vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda+t/2,'t',2,h) *  propag(Lambda,-Lambda,se,dse,'s') * propag(Lambda,-Lambda+t,se,dse,p2)) ;
                    singsc += real(-1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda+t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda+t/2,'t',2,h) *  propag(Lambda,Lambda,se,dse,'s') * propag(Lambda,Lambda+t,se,dse,p2)) ;
                    p1_new = 'e';
                    p2_new = p2;
                    dse_cutoff_up= ffreqs[nw-1]+t/2;dse_cutoff_low= ffreqs[0]+t/2;}

                else if(p2=='k'){

                    singsc = real( -1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,-Lambda-t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,-Lambda-t/2,'t',2,h) *  propag(Lambda,-Lambda-t,se,dse,p1) * propag(Lambda,-Lambda,se,dse,'s')) ;
                    singsc += real(-1./(2*pi)*vert1.vvalsmooth(red_side,map1,a,b,c,t,Lambda-t/2,w2,'t',1,h) * vert2.vvalsmooth(red_side,map2,d,e,f,t,w1,Lambda-t/2,'t',2,h) *  propag(Lambda,Lambda-t,se,dse,p1) * propag(Lambda,Lambda,se,dse,'s')) ;
                    p1_new = p1;
                    p2_new = 'e';
                    dse_cutoff_up= ffreqs[nw-1]-t/2;dse_cutoff_low= ffreqs[0]-t/2;};
                double lhs=0, rhs=0;//result on left and right side of sharp cutoff

                struct tbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension

                if(lower <= limit_low && abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){//left side of sharp cutoff

                    double up_eff=0;
                    if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                    F.params = &params_mod;
                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,dse_cutoff_low ,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfl_re, &errorgfl_re);

                    double compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    lhs=compl_real;

                }
                else if(lower <= limit_up ){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real = (resultgfl_re )*vert1_const * vert2_const;
                    };
                    if(dse_cutoff_low <lower){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                        F.params = &params_mod;
                        F.function = &tbubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                    };

                    lhs=result_re + compl_real;
                }

                else if(lower > limit_up){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <limit_low){
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,dse_cutoff_low,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real = resultgfl_re*vert1_const * vert2_const;

                    };


                    if(dse_cutoff_low <limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &tbubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                    };


                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_low <lower){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                        double up_eff=0;
                        if(dse_cutoff_up < lower){up_eff = dse_cutoff_up;}else{up_eff=lower;};
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real += (resultgfl_re  )*vert1_const * vert2_const;

                    };

                    lhs=result_re + compl_real;
                };

                if(upper >= limit_up && abs(vert1_const * vert2_const)>0 && dse_cutoff_up >upper){//right side of sharp cutoff

                    double low_eff=0;
                    if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                    F.params = &params_mod;
                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                        w, &resultgfu_re, &errorgfu_re);

                    double compl_real = (resultgfu_re  )*vert1_const * vert2_const;

                    rhs=compl_real;
                }
                else if(upper >= limit_low ){
                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up >limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff ,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfu_re  )*vert1_const * vert2_const;

                    };
                    if(dse_cutoff_up > upper){
                        double low_eff=0;
                        if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &tbubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                    };

                    rhs=result_re + compl_real;

                }

                else if(upper < limit_low){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > upper){
                        double low_eff=0;
                        if(dse_cutoff_low > upper){low_eff = dse_cutoff_low;}else{low_eff=upper;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_low){up_eff = dse_cutoff_up;}else{up_eff=limit_low;};
                        F.params = &params_mod;

                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };


                    if(dse_cutoff_up >limit_low){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_low){low_eff = dse_cutoff_low;}else{low_eff=limit_low;};
                        double up_eff=0;
                        if(dse_cutoff_up < limit_up){up_eff = dse_cutoff_up;}else{up_eff=limit_up;};
                        F.params = &params_mod;
                        F.function = &tbubble_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,up_eff,abs_error, rel_error,1500,2,
                                            w, &result_re, &error_re);
                    };


                    if(abs(vert1_const * vert2_const)>0 && dse_cutoff_up > limit_up){
                        double low_eff=0;
                        if(dse_cutoff_low > limit_up){low_eff = dse_cutoff_low;}else{low_eff=limit_up;};

                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,low_eff,dse_cutoff_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real +=(resultgfu_re  )*vert1_const * vert2_const;
                    };

                    rhs=result_re + compl_real;
                };

                B= rhs + lhs + singsc;
            }

            else if(REG==2 && p1=='g' && p2=='g'){

                F.function = &tbubble_re<T1,T2>;
                F.params = &params;



                gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                    w, &result_re, &error_re);
                //integration of greens function outside of w-dependent interval

                double compl_real=0;
                if(abs(vert1_const * vert2_const)>0){
                    F.function = tgreensfunc_re<T1,T2>;
                    gsl_integration_qagil(&F,limit_low,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfl_re, &errorgfl_re);

                    F.function = &tgreensfunc_re<T1,T2>;
                    gsl_integration_qagiu(&F,limit_up,abs_error_bare,rel_error_bare,1500,
                                          w, &resultgfu_re, &errorgfu_re);

                    compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                };
                B=result_re + compl_real ;
            }
            else if(REG==2 && ((p1=='g' && p2=='s') || (p1=='s' && p2=='g'))){

                double bound_low=0, bound_up=0;

                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;

                if(p1=='s'){bound_low = -7*Lambda+t/2; bound_up = 7*Lambda+t/2;}
                else if(p2=='s'){bound_low = -7*Lambda-t/2; bound_up = 7*Lambda-t/2;}

                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;

                    };




                    B=result_re + compl_real ;
                }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    };


                    B=(result_re + compl_real );
                }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){

                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real =(resultgfu_re )*vert1_const * vert2_const;
                    };



                    B=result_re + compl_real ;
                }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real =(resultgfl_re )*vert1_const * vert2_const;
                    };


                    B=( compl_real );
                }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    B=result_re  ;
                }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare,rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };

                    B=( compl_real );
                }

            }

            else if(REG==2 && ((p1=='g' && p2=='k') || (p1=='k' && p2=='g'))){

                double singsc=0;
                double bound_low=0, bound_up=0;
                //FIRST COMPUTE THE SINGLE SCALE CONTRIBUTION
                int m1=0;//turns on and off the integration in the three parts
                int m2=0;
                int m3=0;
                char p1_new,p2_new;
                if(p1=='k'){
                    bound_low = -7*Lambda+t/2; bound_up = 7*Lambda+t/2;
                    p1_new = 's'; p2_new = 'g';}
                else if(p2=='k'){bound_low = -7*Lambda-t/2; bound_up = 7*Lambda-t/2;
                    p1_new = 'g'; p2_new = 's';
                }


                struct tbubble_params<T1,T2> params_singsc= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real =(resultgfl_re +resultgfu_re )*vert1_const * vert2_const;
                    };




                    singsc=result_re + compl_real ;
                }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error, rel_error,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real =(resultgfl_re  )*vert1_const * vert2_const;

                    };


                    singsc=result_re + compl_real ;
                }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){

                        F.params = &params_singsc;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };



                    singsc=result_re + compl_real ;
                }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real =(resultgfl_re )*vert1_const * vert2_const;
                    };


                    singsc=compl_real ;
                }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_singsc;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    singsc=result_re ;
                }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_singsc;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real =(resultgfu_re )*vert1_const * vert2_const;
                    };

                    singsc = compl_real ;
                }
                //ADD THE CONTRIBUTION FROM THE KATANON EXTENSION:
                double kat=0;

                //reset m-values
                m1=0;//turns on and off the integration in the three parts
                m2=0;
                m3=0;


                if(p1=='k'){
                    bound_low = ffreqs[0]+t/2; bound_up = ffreqs[nw-1]+t/2;//the bounds are the boundaries of the dse-grid
                    p1_new = 'e';
                    p2_new = 'g';}
                else if(p2=='k'){
                    bound_low = ffreqs[0]-t/2; bound_up = ffreqs[nw-1]-t/2;
                    p1_new = 'g';
                    p2_new = 'e';}

                struct tbubble_params<T1,T2> params_mod= {red_side,map1,map2,Lambda,vert1, vert2,p1_new,p2_new,se, dse, a,  b,  c,  d,  e,  f,  t, w1,  w2,h};//modify integration parameters such that one of the propagators only yields the katanin extension


                if(bound_low < limit_low){
                    m1=1;
                    if(bound_up > limit_low){
                        m2=1;
                        if(bound_up > limit_up){
                            m3=1;
                        };
                    };
                }
                else if(bound_low < limit_up){
                    m2=1;
                    if(bound_up > limit_up){
                        m3=1;
                    };
                }
                else{m3=1;};

                if(m1==1 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfl_re +resultgfu_re)*vert1_const * vert2_const;

                    };




                    kat=result_re + compl_real ;
                }


                else if(m1==1 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,limit_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,limit_low,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);
                        compl_real = (resultgfl_re  )*vert1_const * vert2_const;

                    };


                    kat=result_re + compl_real ;
                }

                else if(m1==0 && m2==1 && m3==1){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,limit_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);
                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){


                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,limit_up,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);

                        compl_real = (resultgfu_re )*vert1_const * vert2_const;

                    };



                    kat=result_re + compl_real ;
                }

                else if(m1==1 && m2==0 && m3==0){


                    //integration of greens function outside of w-dependent interval

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfl_re, &errorgfl_re);

                        compl_real = (resultgfl_re)*vert1_const * vert2_const;

                    };


                    kat= compl_real ;
                }
                else if(m1==0 && m2==1 && m3==0){

                    F.function = &tbubble_re<T1,T2>;
                    F.params = &params_mod;

                    gsl_integration_qag(&F,bound_low,bound_up,abs_error, rel_error,1500,2,
                                        w, &result_re, &error_re);

                    //integration of greens function outside of w-dependent interval

                    kat=result_re  ;
                }
                else if(m1==0 && m2==0 && m3==1){

                    double compl_real=0;
                    if(abs(vert1_const * vert2_const)>0){
                        F.params = &params_mod;
                        F.function = &tgreensfunc_re<T1,T2>;
                        gsl_integration_qag(&F,bound_low,bound_up,abs_error_bare, rel_error_bare,1500,2,
                                            w, &resultgfu_re, &errorgfu_re);
                        compl_real =(resultgfu_re )*vert1_const * vert2_const;

                    };

                    kat= compl_real ;
                }
                B=kat + singsc;

            };


        };


    };



    if(abs(B)<1e-20){B=0;};//prevents errors from unprecise cancellation of diagrams
    return B;
};






#endif //KELDYSH_MFRG_BUBBLES_H

#pragma clang diagnostic pop