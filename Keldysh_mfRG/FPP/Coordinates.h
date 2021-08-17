//
// Created by Marcel on 14.04.2021.
//

#ifndef MAIN_CPP_COORDINATES_H
#define MAIN_CPP_COORDINATES_H

#include <iostream>
#include <complex>          // for usage of complex numbers
#include <cmath>            // for math. operations (real, imag, abs etc.)
#include <vector>           // vec class is derived from vector class
#include <initializer_list> // to initialize vec class with initializer list


std::vector<double> cartesion_to_spherical(std::vector<double> v){
    double x = v[0];
    double y = v[1];
    double z = v[2];
    double r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    double theta;
    if (z != 0) {
        theta = atan(sqrt(pow(x, 2) + pow(y, 2)) / z);
    }    else {
        theta = M_PI / 2;
    }
    if (theta < 0) {
        theta = theta + 2*M_PI;
    }
    double phi;
    if (x > 0) {
        phi = atan(y / x);
    }   else if(x < 0 && y >= 0) {
        phi = atan(y / x) + M_PI;
    }
    else if(x < 0 && y < 0) {
        phi = atan(y / x) - M_PI;
    }
    else if(x == 0 && y > 0) {
        phi = M_PI / 2;
    }   else if(x == 0 && y < 0) {
        phi = -M_PI / 2;
    }   else {
        phi = 0;
    }
    if (phi < 0) {
        phi = phi + 2*M_PI;
    }
    std::vector<double> vsph;
    vsph = {r, theta, phi};
    return vsph;
}

std::vector<double> spherical_to_cartesion(std::vector<double> v) {
    double r = v[0];
    double theta = v[1];
    double phi = v[2];
    double x = r * sin(theta) * cos(phi);
    double y = r * sin(theta) * sin(phi);
    double z = r * cos(theta);
    std::vector<double> vcar;
    vcar = {x, y, z};
    return vcar;
}

std::vector<double> transform_back_forth(std::vector<double> v, int N) {
    std::vector<double> v_old;
    std::vector<double> v_new;
    v_old = v;
    for (int i = 0; i < N; ++i){
        v_new = cartesion_to_spherical(v_old);
        v_old = spherical_to_cartesion(v_new);
    }
    return v_old;
}

double test_for_loop(double x, int N) {
    double x_old;
    double x_new;
    x_old = x;
    for (int i = 0; i < N; ++i){
        x_new = x_old + 2.0;
        x_old = x_new + 2.0;
    }
    return x_old;
}

void coordinate_transform(){
    int transform;
    std::cout << "If you want to transform from cartesian to spherical coordinates, type in '0', if you want to transform from spherical coordinates to cartesian, type in '1'.\n";
    std::cin >> transform;
    std::vector<double> v_in;
    std::vector<double> v_out;
    double v1_in;
    double v2_in;
    double v3_in;
    double v1_out;
    double v2_out;
    double v3_out;
    if (transform == 0) {
        std::cout << "Type in x-component of a vector in cartesian coordinates: \n";
        std::cin >> v1_in;
        std::cout << "Type in y-component of a vector in cartesian coordinates: \n";
        std::cin >> v2_in;
        std::cout << "Type in z-component of a vector in cartesian coordinates: \n";
        std::cin >> v3_in;
        v_in = {v1_in, v2_in, v3_in};
        v_out = cartesion_to_spherical(v_in);
        v1_out = v_out[0];
        v2_out = v_out[1];
        v3_out = v_out[2];
        std::cout << "In spherical coordinates your vector is {" << v1_out << "," << v2_out << "," << v3_out << "}.\n";
    } else if (transform == 1) {
        std::cout << "Type in radial component of a vector in spherical coordinates: \n";
        std::cin >> v1_in;
        std::cout << "Type in polar angle of a vector in spherical coordinates: \n";
        std::cin >> v2_in;
        std::cout << "Type in azimuthal angle of a vector in spherical coordinates: \n";
        std::cin >> v3_in;
        v_in = {v1_in, v2_in, v3_in};
        v_out = spherical_to_cartesion(v_in);
        v1_out = v_out[0];
        v2_out = v_out[1];
        v3_out = v_out[2];
        std::cout << "In cartesian coordinates your vector is {" << v1_out << "," << v2_out << "," << v3_out << "}.\n";
    } else {
        int n;
        std::cout << "Type in an integer number for back and forth iterations.\n";
        std::cin >> n;
        std::cout << "Type in x-component of a vector in cartesian coordinates: \n";
        std::cin >> v1_in;
        std::cout << "Type in y-component of a vector in cartesian coordinates: \n";
        std::cin >> v2_in;
        std::cout << "Type in z-component of a vector in cartesian coordinates: \n";
        std::cin >> v3_in;
        v_in = {v1_in, v2_in, v3_in};
        v_out = transform_back_forth(v_in, n);
        v1_out = v_out[0];
        v2_out = v_out[1];
        v3_out = v_out[2];
        std::cout << "After " << n << " back and forth transformations your vector is {" << v1_out << "," << v2_out << "," << v3_out << "}.\n";
    }
}

#endif //MAIN_CPP_COORDINATES_H
