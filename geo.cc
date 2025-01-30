#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TLine.h>
#include <tuple>
#include <TGraph.h>
#include <vector>
#include <cmath>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <stdexcept>
#include "TVectorD.h"
#include <fstream>
#include <random>
#include <thread>
#include <chrono>
#include<TF1.h>

#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewCell.hh"

using namespace Garfield;
using namespace std;
using Point3D = array<double, 3>;
using Matrix3x3 = array<array<double, 3>, 3>;
const double c = 299792458.0;

const double rCell = 0.5;
double offset = 6.35;
double length = 11.5;
double girth = 17.2;
double sig = 0;
int layers = 4;
int rows = static_cast<int>(length / (2*rCell));

TH2F* pos_z = new TH2F("pos_z","pos_z",30,-7.5,7.5,10,-2.5,2.5);
TH1F* reco_z = new TH1F("reco_z","reco_z",1000,-25,25);
TH2F* reco_z_hits = new TH2F("reco_z_hits","reco_z_hits",10,3,13,1000,-25,25);

struct Cylinder {
    double a, b, c;  // Center position (a, b, c)
    double r;        // Radius
    double L;        // Length
    int mod;         // Module

    // Constructor
    Cylinder(double x, double y, double z, double radius, double length, int stack)
        : a(x), b(y), c(z), r(radius), L(length), mod(stack) {}
};

double randomiser(double radius) {
    static std::random_device rd; // Random device for seeding
    static std::mt19937 gen(rd()); // Mersenne Twister generator
    static std::uniform_real_distribution<> dis(0.9,1.1);

    double rad = dis(gen);

    if (rad*radius>rCell) {return rCell;}
    else {return rad*radius;}
}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCentersLongitudinal() {

    std::vector<double> xCenter, yCenter, zCenter;

    double offset = 2.5;

    double Vspace = 2.5;
    double Hspace = 10;

    int nTubesV = static_cast<int>(0.5*Vspace / rCell);
    int nTubesH = static_cast<int>(0.5*Hspace / rCell);

    //Original
    for (int j = 0; j<nTubesH; ++j) {
        for (int i = -j-1; i <= 2*nTubesV+j; ++i) {
                xCenter.push_back(offset+j*2*rCell);
                yCenter.push_back(-nTubesV*2*rCell+i*2*rCell+rCell);
                zCenter.push_back(0);
        }
    }

    //Rotate 90
    std::vector<double> x_rotated,y_rotated,z_rotated;

    for (int j = 0; j<nTubesH; ++j) {
        for (int i = -j; i <= 2*nTubesV+j-1; ++i) {
                y_rotated.push_back(-(offset+j*2*rCell));
                x_rotated.push_back(-nTubesV*2*rCell+i*2*rCell+rCell);
                z_rotated.push_back(0);
        }
    }

    //Reflect verticals
    std::vector<double> x_flip,y_flip,z_flip;

    for (size_t i = 0; i < xCenter.size(); ++i) {
        x_flip.push_back(-xCenter[i]);
        y_flip.push_back(yCenter[i]);
        z_flip.push_back(zCenter[i]);  // z remains unchanged
    }

      // Append the rotated coordinates to the original vectors
    xCenter.insert(xCenter.end(), x_flip.begin(), x_flip.end());
    yCenter.insert(yCenter.end(), y_flip.begin(), y_flip.end());
    zCenter.insert(zCenter.end(), z_flip.begin(), z_flip.end());

    //Add rotated
    for (size_t i = 0; i < x_rotated.size(); ++i) {
        xCenter.push_back(x_rotated[i]);
        yCenter.push_back(y_rotated[i]);
        zCenter.push_back(z_rotated[i]);  // z remains unchanged
    }


    for (size_t i = 0; i < x_rotated.size(); ++i) {
        xCenter.push_back(x_rotated[i]);
        yCenter.push_back(-y_rotated[i]);
        zCenter.push_back(z_rotated[i]);  // z remains unchanged
    }

    return std::make_tuple(xCenter, yCenter, zCenter);

}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCentersTransverseStep() {

    std::vector<double> xCenter, yCenter, zCenter;

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-offset-1);
            yCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
            zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
            yCenter.push_back(offset+1);
            zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset+1);
            yCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
            zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);
        }
    }


    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
            yCenter.push_back(-offset-1);
            zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);
        }
    }

    return std::make_tuple(xCenter, yCenter, zCenter);

}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>,std::vector<double>> GetStrawCentersTransverseStepZoffsetOverhang() {

    std::vector<double> xCenter, yCenter, zCenter;
    std::vector<double> mod_straws;
    double count = 0;

    for (int j = 0; j<layers; ++j) {

        if (j==0) {
            for (int i = 0; i<21; i++) {
                xCenter.push_back(offset-girth/2);
                yCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
                zCenter.push_back(13-0.1-i*0.1-rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==1) {
            for (int i = 0; i<20; i++) {
                xCenter.push_back(offset-girth/2);
                yCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
                zCenter.push_back(13-0.15-i*0.1-2*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==2) {
            for (int i = 0; i<19; i++) {
                xCenter.push_back(offset-girth/2);
                yCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
                zCenter.push_back(13-0.25-i*0.1-3*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==3) {
            for (int i = 0; i<17; i++) {
                xCenter.push_back(offset-girth/2);
                yCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
                zCenter.push_back(13-0.40-i*0.1-6*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

    }

    cout<<count<<endl;
    mod_straws.push_back(count);

    for (int j = 0; j<layers; ++j) {

        if (j==0) {
            for (int i = 0; i<21; i++) {
                xCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
                yCenter.push_back(-offset+girth/2);
                zCenter.push_back(13-0.1-i*0.1-rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==1) {
            for (int i = 0; i<20; i++) {
                xCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
                yCenter.push_back(-offset+girth/2);
                zCenter.push_back(13-0.15-i*0.1-2*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==2) {
            for (int i = 0; i<19; i++) {
                xCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
                yCenter.push_back(-offset+girth/2);
                zCenter.push_back(13-0.25-i*0.1-3*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==3) {
            for (int i = 0; i<17; i++) {
                xCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
                yCenter.push_back(-offset+girth/2);
                zCenter.push_back(13-0.4-i*0.1-6*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

    }

    cout<<count<<endl;
    mod_straws.push_back(count);

    for (int j = 0; j<layers; ++j) {

        if (j==0) {
            for (int i = 0; i<21; i++) {
                xCenter.push_back(girth/2-offset);
                yCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
                zCenter.push_back(13-0.1-i*0.1-rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==1) {
            for (int i = 0; i<20; i++) {
                xCenter.push_back(girth/2-offset);
                yCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
                zCenter.push_back(13-0.15-i*0.1-2*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==2) {
            for (int i = 0; i<19; i++) {
                xCenter.push_back(girth/2-offset);
                yCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
                zCenter.push_back(13-0.25-i*0.1-3*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==3) {
            for (int i = 0; i<17; i++) {
                xCenter.push_back(girth/2-offset);
                yCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
                zCenter.push_back(13-0.4-i*0.1-6*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

    }

    cout<<count<<endl;
    mod_straws.push_back(count);

    for (int j = 0; j<layers; ++j) {

        if (j==0) {
            for (int i = 0; i<21; i++) {
                xCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
                yCenter.push_back(+offset-girth/2);
                zCenter.push_back(13-0.1-i*0.1-rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==1) {
            for (int i = 0; i<20; i++) {
                xCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
                yCenter.push_back(+offset-girth/2);
                zCenter.push_back(13-0.15-i*0.1-2*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==2) {
            for (int i = 0; i<19; i++) {
                xCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
                yCenter.push_back(+offset-girth/2);
                zCenter.push_back(13-0.25-i*0.1-3*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

        if (j==3) {
            for (int i = 0; i<17; i++) {
                xCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
                yCenter.push_back(+offset-girth/2);
                zCenter.push_back(13-0.4-i*0.1-6*rCell-i*2*rCell);
                count=count+1;
                //cout<<"("<<xCenter.back()<<","<<yCenter.back()<<","<<zCenter.back()<<") "<<count<<endl;
            }
        }

    }

    cout<<count<<endl;
    mod_straws.push_back(count);
    cout<<mod_straws[0]<<","<<mod_straws[1]<<","<<mod_straws[2]<<","<<mod_straws[3]<<endl;

    return std::make_tuple(xCenter, yCenter, zCenter,mod_straws);

}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCentersTransverseStepZoffset() {

    std::vector<double> xCenter, yCenter, zCenter;
    int count = 0;
    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset-girth/2);
            yCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
            //zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);

            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell+i*0.1+0.1);}
            else {zCenter.push_back(i*2*rCell+i*0.1+0.05);}
            count = count + 1;
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
            yCenter.push_back(-offset+girth/2);
            //zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);

            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell+i*0.1+0.1);}
            else {zCenter.push_back(i*2*rCell+i*0.1+0.05);}
            count = count + 1;
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(girth/2-offset);
            yCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
            //zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);

            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell+i*0.1+0.1);}
            else {zCenter.push_back(i*2*rCell+i*0.1+0.05);}
            count = count + 1;
        }
    }


    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
            yCenter.push_back(+offset-girth/2);
            //zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);

            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell+i*0.1+0.1);}
            else {zCenter.push_back(i*2*rCell+i*0.1+0.05);}
            count = count + 1;
        }
    }

    cout<<count<<endl;

    return std::make_tuple(xCenter, yCenter, zCenter);

}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCentersTransverseStepXYoffset() {

    std::vector<double> xCenter, yCenter, zCenter;

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            //xCenter.push_back(-offset-1);
            yCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
            zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);

            if (j%2!=0) {xCenter.push_back(-offset-rCell-1);}
            else {xCenter.push_back(-offset-1);}
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
            //yCenter.push_back(offset+1);
            zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);

            if (j%2!=0) {yCenter.push_back(offset+1);}
            else {yCenter.push_back(offset+rCell+1);}
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            //xCenter.push_back(offset+1);
            yCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
            zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);


            if (j%2!=0) {xCenter.push_back(offset+rCell+1);}
            else {xCenter.push_back(offset+1);}
        }
    }


    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
            //yCenter.push_back(-offset-1);
            zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);

            if (j%2!=0) {yCenter.push_back(-offset-1);}
            else {yCenter.push_back(-offset-rCell-1);}
        }
    }

    return std::make_tuple(xCenter, yCenter, zCenter);

}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCentersTransverseStepZXYoffset() {

    std::vector<double> xCenter, yCenter, zCenter;

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            //xCenter.push_back(-offset-1);
            yCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
            //zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);

            if (j%2!=0) {xCenter.push_back(-offset-rCell-1);}
            else {xCenter.push_back(-offset-1);}

            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);}
            else {zCenter.push_back(i*2*rCell+i*0.2+0.1);}
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset+j*(2*rCell+0.1)+rCell+0.1);
            //yCenter.push_back(offset+1);
            //zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);

            if (j%2!=0) {yCenter.push_back(offset+1);}
            else {yCenter.push_back(offset+rCell+1);}

            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);}
            else {zCenter.push_back(i*2*rCell+i*0.2+0.1);}
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            //xCenter.push_back(offset+1);
            yCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
            //zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);


            if (j%2!=0) {xCenter.push_back(offset+rCell+1);}
            else {xCenter.push_back(offset+1);}

            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);}
            else {zCenter.push_back(i*2*rCell+i*0.2+0.1);}
        }
    }


    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-(offset+j*(2*rCell+0.1)+rCell+0.1));
            //yCenter.push_back(-offset-1);
            //zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);

            if (j%2!=0) {yCenter.push_back(-offset-1);}
            else {yCenter.push_back(-offset-rCell-1);}

            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell+i*0.2+0.2);}
            else {zCenter.push_back(i*2*rCell+i*0.2+0.1);}
        }
    }

    return std::make_tuple(xCenter, yCenter, zCenter);

}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCentersTransverse() {

    std::vector<double> xCenter, yCenter, zCenter;

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-offset);
            yCenter.push_back(offset+j*2*rCell+rCell);
            zCenter.push_back(i*2*rCell+rCell);
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset+j*2*rCell+rCell);
            yCenter.push_back(offset);
            zCenter.push_back(i*2*rCell+rCell);
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset);
            yCenter.push_back(-(offset+j*2*rCell+rCell));
            zCenter.push_back(i*2*rCell+rCell);
        }
    }



    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-(offset+j*2*rCell+rCell));
            yCenter.push_back(-offset);
            zCenter.push_back(i*2*rCell+rCell);
        }
    }



    return std::make_tuple(xCenter, yCenter, zCenter);

}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCentersTransverseSingleOffset() {

     std::vector<double> xCenter, yCenter, zCenter;

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-offset);
            yCenter.push_back(offset+j*2*rCell+rCell);
            //zCenter.push_back(i*2*rCell+rCell);
            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            else {zCenter.push_back(i*2*rCell);}
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset+j*2*rCell+rCell);
            yCenter.push_back(offset);
            //zCenter.push_back(i*2*rCell+rCell);
            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            else {zCenter.push_back(i*2*rCell);}
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset);
            yCenter.push_back(-(offset+j*2*rCell+rCell));
            //zCenter.push_back(i*2*rCell+rCell);
            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            else {zCenter.push_back(i*2*rCell);}
        }
    }



    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-(offset+j*2*rCell+rCell));
            yCenter.push_back(-offset);
            //zCenter.push_back(i*2*rCell+rCell);
            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            else {zCenter.push_back(i*2*rCell);}
        }
    }



    return std::make_tuple(xCenter, yCenter, zCenter);


}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCentersTransverseDoubleOffset() {

    std::vector<double> xCenter, yCenter, zCenter;

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            //xCenter.push_back(-offset);
            yCenter.push_back(offset+j*2*rCell+rCell);
            //zCenter.push_back(i*2*rCell+rCell);
            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            else {zCenter.push_back(i*2*rCell);}

            if (j%2!=0) {xCenter.push_back(-offset-rCell);}
            else {xCenter.push_back(-offset);}
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset+j*2*rCell+rCell);
            //yCenter.push_back(offset);
            //zCenter.push_back(i*2*rCell+rCell);
            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            else {zCenter.push_back(i*2*rCell);}

            if (j%2!=0) {yCenter.push_back(offset);}
            else {yCenter.push_back(offset+rCell);}
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            //xCenter.push_back(offset);
            yCenter.push_back(-(offset+j*2*rCell+rCell));
            //zCenter.push_back(i*2*rCell+rCell);
            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            else {zCenter.push_back(i*2*rCell);}

            if (j%2!=0) {xCenter.push_back(offset+rCell);}
            else {xCenter.push_back(offset);}
        }
    }



    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-(offset+j*2*rCell+rCell));
            yCenter.push_back(-offset);
            //zCenter.push_back(i*2*rCell+rCell);
            if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            else {zCenter.push_back(i*2*rCell);}

            if (j%2!=0) {yCenter.push_back(-offset);}
            else {yCenter.push_back(-offset-rCell);}
        }
    }



    return std::make_tuple(xCenter, yCenter, zCenter);

}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCentersTransverseXYOffset() {

    std::vector<double> xCenter, yCenter, zCenter;

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            //xCenter.push_back(-offset);
            yCenter.push_back(offset+j*2*rCell+rCell);
            zCenter.push_back(i*2*rCell+rCell);
            //if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            //else {zCenter.push_back(i*2*rCell);}

            if (j%2!=0) {xCenter.push_back(-offset-rCell);}
            else {xCenter.push_back(-offset);}
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset+j*2*rCell+rCell);
            //yCenter.push_back(offset);
            zCenter.push_back(i*2*rCell+rCell);
            //if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            //else {zCenter.push_back(i*2*rCell);}

            if (j%2!=0) {yCenter.push_back(offset);}
            else {yCenter.push_back(offset+rCell);}
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            //xCenter.push_back(offset);
            yCenter.push_back(-(offset+j*2*rCell+rCell));
            zCenter.push_back(i*2*rCell+rCell);
            //if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            //else {zCenter.push_back(i*2*rCell);}

            if (j%2!=0) {xCenter.push_back(offset+rCell);}
            else {xCenter.push_back(offset);}
        }
    }



    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-(offset+j*2*rCell+rCell));
            yCenter.push_back(-offset);
            //zCenter.push_back(i*2*rCell+rCell);
            //if (j%2==0) {zCenter.push_back(i*2*rCell+rCell);}
            //else {zCenter.push_back(i*2*rCell);}

            if (j%2!=0) {yCenter.push_back(-offset);}
            else {yCenter.push_back(-offset-rCell);}
        }
    }



    return std::make_tuple(xCenter, yCenter, zCenter);

}

void PlotMidYZ(const std::vector<double>& zCenter, const std::vector<double>& yCenter, std::vector<double> mod_straws) {
    TCanvas* canvas = new TCanvas("canvas", "Straws and Track", 800, 800);
    // Create a canvas
    canvas->SetFixedAspectRatio(); // Ensure equal scaling on both axes
    canvas->DrawFrame(-20, -20, 20, 20); // Set the drawing frame (xMin, yMin, xMax, yMax)

    TEllipse* straw = new TEllipse(0,0,15);
            straw->SetFillStyle(0); // No fill
            straw->SetLineColor(kBlack);
            straw->Draw("same");

    TBox* cell = new TBox(-7.5, -3, 7.5, 3);
        cell->SetFillStyle(0); // No fill
        cell->SetLineColor(kRed);
        cell->Draw("same");

    TBox* vacuum = new TBox(-13.5, -6.35, 13.5, 6.35);
        vacuum->SetFillStyle(0); // No fill
        vacuum->SetLineColor(kBlue);
        vacuum->Draw("same");

    // Draw the circles (straws)
    for (size_t i = 0; i < yCenter.size(); ++i) {

        if (i<mod_straws[0] || (i>=mod_straws[1] && i<mod_straws[2]) ) {
            TEllipse* straw = new TEllipse(zCenter[i], yCenter[i], rCell);
                straw->SetFillStyle(0); // No fill
                straw->SetLineColor(kBlack);
                straw->Draw("same");
        }

    }

    canvas->Update();
}

void PlotDOCAZY(const std::vector<double>& zcell, const std::vector<double>& ycell, std::vector<double> radius, std::vector<double> mod_straws) {

    // Draw the circles (straws)
    for (size_t i = 0; i < radius.size(); ++i) {

        if (i<mod_straws[0] || (i>mod_straws[1] && i<mod_straws[2]) ) {
            TEllipse* DOCA = new TEllipse(zcell[i], ycell[i], radius[i]);
            DOCA->SetFillStyle(0); // No fill
            DOCA->SetLineColor(kRed);
            DOCA->Draw("same");
        }
    }

}

void PlotXYCells(const std::vector<double>& xCenter, const std::vector<double>& yCenter, std::vector<double> mod_straws) {
    // Check if the inputs are valid

    TCanvas* canvas = new TCanvas("canvas", "Straws and Track", 1000, 1000);
    // Create a canvas
    canvas->SetFixedAspectRatio(); // Ensure equal scaling on both axes
    canvas->DrawFrame(-20, -20, 20, 20); // Set the drawing frame (xMin, yMin, xMax, yMax)

    TEllipse* DOCA = new TEllipse(0,0,15);
            DOCA->SetFillStyle(0); // No fill
            DOCA->SetLineColor(kRed);
            DOCA->Draw("same");

    TEllipse* DOCA1 = new TEllipse(0,0,6.35);
            DOCA1->SetFillStyle(0); // No fill
            DOCA1->SetLineColor(kBlack);
            DOCA1->Draw("same");


    // Draw the box (straws)
    for (int i = 0; i < mod_straws[0]; ++i) {
        TBox* straw = new TBox(xCenter[i]-girth/2+3, yCenter[i]-rCell, xCenter[i]+girth/2-1, yCenter[i]+rCell);
        straw->SetFillStyle(0); // No fill
        straw->SetLineColor(kBlack);
        straw->Draw("same");
    }

    for (int i = mod_straws[0]; i < mod_straws[1]; ++i) {
        TBox* straw1 = new TBox(xCenter[i]-rCell, yCenter[i]-girth/2+1, xCenter[i]+rCell, yCenter[i]+girth/2-3);
        straw1->SetFillStyle(0); // No fill
        straw1->SetLineColor(kRed);
        straw1->Draw("same");
    }

    for (int i = mod_straws[1]; i < mod_straws[2]; ++i) {
        TBox* straw2 = new TBox(xCenter[i]-girth/2+1, yCenter[i]-rCell, xCenter[i]+girth/2-3, yCenter[i]+rCell);
        straw2->SetFillStyle(0); // No fill
        straw2->SetLineColor(kBlue);
        straw2->Draw("same");
    }

    for (int i = mod_straws[2]; i < mod_straws[3]; ++i) {
        TBox* straw3 = new TBox(xCenter[i]-rCell, yCenter[i]-girth/2+3, xCenter[i]+rCell, yCenter[i]+girth/2-1);
        straw3->SetFillStyle(0); // No fill
        straw3->SetLineColor(kGreen);
        straw3->Draw("same");

    }


    // Update the canvas to render the drawing
    canvas->Update();
}

void PlotTrack(std::vector<double> trk_x,std::vector<double> trk_y) {

  TGraph* graph = new TGraph(trk_x.size(), trk_x.data(), trk_y.data());

  // Set the marker style for the points (e.g., a circle)
  graph->SetMarkerStyle(21);  // 21: circle
  graph->SetMarkerColor(kRed);  // Red points
  graph->SetMarkerSize(3);  // Marker size

  // Draw the graph on the canvas
  graph->Draw("same");

}

void PlotFit(std::vector<double> trk_x,std::vector<double> trk_y) {

  TGraph* graph = new TGraph(trk_x.size(), trk_x.data(), trk_y.data());

  // Set the marker style for the points (e.g., a circle)
  graph->SetMarkerStyle(21);  // 21: circle
  graph->SetLineColor(kBlue);  // Red points
  graph->SetMarkerSize(3);  // Marker size

  // Draw the graph on the canvas
  graph->Draw("same");

}

std::vector<double> generateDir() {
    // Seed the random number generator
    srand(static_cast<unsigned>(time(nullptr)));

    static std::random_device rd; // Random device for seeding
    static std::mt19937 gen(rd()); // Mersenne Twister generator
    static std::uniform_real_distribution<> dis(0,2);

    // Generate uniform random angles
    double phi = M_PI*dis(gen);  // Azimuthal angle (0 to 2π)
    double theta = acos(1.0 - dis(gen));  // Polar angle (0 to π)

    // Convert spherical coordinates to Cartesian coordinates
    double x = sin(theta) * cos(phi);
    double y = sin(theta) * sin(phi);
    double z = cos(theta);

    return {x,y,z};

}

std::vector<double> generatePos() {
    // Seed the random number generator
    srand(static_cast<unsigned>(time(nullptr)));

    static std::random_device rd; // Random device for seeding
    static std::mt19937 gen(rd()); // Mersenne Twister generator
    static std::uniform_real_distribution<> disXY(-0.3,0.3);
    static std::uniform_real_distribution<> disZ(-7.5,7.5);

    // Generate uniform random angles
    double x = disXY(gen);
    double y = disXY(gen);
    double z = disZ(gen);

    //cout<<x<<" "<<y<<" "<<z<<endl;

    return {x,y,z};

}

std::vector<double> vertex_gen() {

    // Seed the random number generator
    srand(static_cast<unsigned>(time(nullptr)));

    static std::random_device rd; // Random device for seeding
    static std::mt19937 gen(rd()); // Mersenne Twister generator
    static std::uniform_real_distribution<> disXY(-0.3,0.3);//XY dispursion
    static std::uniform_real_distribution<> disAB(-5.0/1000.0,5.0/1000.0);//Angular dispursion
    static std::uniform_real_distribution<> disZ(-2.5,2.5);//Z vertex

    // Generate incoming angle and position
    double x0 = disXY(gen), y0 = disXY(gen);
    double A = disAB(gen), B = disAB(gen);
    double z0 = sqrt(1.0 - tan(A)*tan(A) - tan(B)*tan(B));
    std::vector<double> direction = {tan(A),tan(B),z0};

    //Calculate incoming direction and positions
    double z = disZ(gen);
    double x = x0 + z*tan(A), y = y0 + z*tan(B);
    std::vector<double> position = {x,y,z};

    return position;

}


double lorentzFactor(double v) {
    return 1.0 / std::sqrt(1.0 - (v * v));
}

// Function to compute the Lorentz boost for direction vector along the z-axis
std::vector<double> lorentzBoostDirection(const std::vector<double>& particleDirection,
                                             double sourceVelocity) {
    // Calculate the Lorentz factor
    double gamma = lorentzFactor(sourceVelocity);

    // Direction components of the particle in the source frame
    double v_x = particleDirection[0];
    double v_y = particleDirection[1];
    double v_z = particleDirection[2];

    // Boosted z-component using relativistic velocity addition
    double v_z_prime = (v_z + sourceVelocity) / (1 + (v_z * sourceVelocity));
    double v_x_prime = v_x;
    double v_y_prime = v_y;

    // Normalize the direction vector
    double magnitude = std::sqrt(v_x_prime * v_x_prime + v_y_prime * v_y_prime + v_z_prime * v_z_prime);
    v_x_prime /= magnitude;
    v_y_prime /= magnitude;
    v_z_prime /= magnitude;

    return {v_x_prime, v_y_prime, v_z_prime};
}


std::pair<std::vector<double>,std::vector<double>> scatter() {

    std::vector<double> position = vertex_gen();

    float knock_m = 938;
    float knock_k = 200;
    float knock_E = knock_k+knock_m;
    float knock_lorenz = knock_E/knock_m;
    float knock_v_mag = sqrt(1-1/(knock_lorenz*knock_lorenz));

    //COM isotropic emmission
    std::vector<double> COM_dir = generateDir();
    std::vector<double> lab_dir = lorentzBoostDirection(COM_dir, knock_v_mag);

    //cout<<"("<<COM_dir[0]<<","<<COM_dir[1]<<","<<COM_dir[2]<<")"<<endl;
    //cout<<"("<<lab_dir[0]<<","<<lab_dir[1]<<","<<lab_dir[2]<<")"<<endl;

    return std::make_pair(lab_dir, position);

}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>,
           std::vector<double>, std::vector<double>, std::vector<double>, std::vector<int>>
hits(std::vector<double> posvect, std::vector<double> dirvect, std::vector<double> xcell, std::vector<double> ycell, std::vector<double> zcell, std::vector<double> mod_straws) {

    std::vector<double> xcells, ycells, zcells;
    std::vector<double> xhits, yhits, zhits;
    std::vector<double> radius;
    std::vector<int> axis;

    double mag = sqrt(dirvect[0]*dirvect[0]+dirvect[1]*dirvect[1]+dirvect[2]*dirvect[2]);
    double ux = dirvect[0]/mag, uy = dirvect[1]/mag, uz = dirvect[2]/mag;
    double A=0,B=0,C=0;
    double outer = 0;
    double inner = 0;

    for (size_t i = 0; i < zcell.size(); i++) {

        double diff1,diff2,x,y,z,t,x1,y1,z1,x2,y2,z2;

        double dx = posvect[0]-xcell[i];
        double dy = posvect[1]-ycell[i];
        double dz = posvect[2]-zcell[i];

        if ((i < mod_straws[0]) || (( (i > mod_straws[1]-1) && (i < mod_straws[2])))) {
            A = uz*uz+uy*uy;
            B = 2*(dz*uz+dy*uy);
            C = dz*dz + dy*dy - rCell*rCell;
        }

        if (((i > mod_straws[0]-1) && (i < mod_straws[1])) || ((i > mod_straws[2]-1) && (i < mod_straws[3]))) {
            A = ux*ux+uz*uz;
            B = 2*(dx*ux+dz*uz);
            C = dx*dx + dz*dz - rCell*rCell;
        }



        if (B*B-4*A*C>0) {

            double t1 = (-B + sqrt(B*B-4*A*C)) / (2 * A);
            double t2 = (-B - sqrt(B*B-4*A*C)) / (2 * A);


            if (t1>0) {
                x1 = posvect[0]+ t1*ux;
                y1 = posvect[1]+ t1*uy;
                z1 = posvect[2]+ t1*uz;

                diff1 = sqrt(x1*x1+y1*y1+z1*z1);

            }

            if (t2>0) {
                x2 = posvect[0]+ t2*ux;
                y2 = posvect[1]+ t2*uy;
                z2 = posvect[2]+ t2*uz;

                diff2 = sqrt(x2*x2+y2*y2+z2*z2);
            }

            if ((diff2>diff1 && diff1!=0)) {x=x1;y=y1;z=z1;}
            if ((diff2<diff1 && diff2!=0)) {x=x2;y=y2;z=z2;}

            double midx=(x1+x2)/2,midy=(y1+y2)/2,midz=(z1+z2)/2;

            //cout<<t1<<" "<<t2<<" "<<x<<" "<<y<<" "<<z<<" "<<A<<" "<<B<<" "<<C<<" "<<B*B-4*A*C<<endl;


            // Draw the box (straws)
    if (i < mod_straws[0]) {
        if (x>= xcell[i] - 2*offset - outer && x<= xcell[i] + 2*offset + inner) {

            if (t1>0 && t2>0) {

                xhits.push_back(x);
                xcells.push_back(xcell[i]);

                yhits.push_back(y);
                ycells.push_back(ycell[i]);

                zhits.push_back(z);
                zcells.push_back(zcell[i]);

                radius.push_back(sqrt((midy-ycell[i])*(midy-ycell[i])+(midz-zcell[i])*(midz-zcell[i])));

                axis.push_back(0);

                }

                //cout<<x<<" "<<y<<" "<<z<<" "<<t1<<" "<<t2<<" "<<xcell[i]<<" "<<ycell[i]<<" "<<zcell[i]<<" "<<sqrt((y-ycell[i])*(y-ycell[i])+(z-zcell[i])*(z-zcell[i]))<<" "<<i<<" codex"<<endl;

            }
        }

    if ( (i > mod_straws[1]-1) && (i < mod_straws[2])) {
        if (x>= xcell[i] - 2*offset - inner && x<= xcell[i] + 2*offset + outer) {

            if (t1>0 && t2>0) {

                xhits.push_back(x);
                xcells.push_back(xcell[i]);

                yhits.push_back(y);
                ycells.push_back(ycell[i]);

                zhits.push_back(z);
                zcells.push_back(zcell[i]);

                radius.push_back(sqrt((midy-ycell[i])*(midy-ycell[i])+(midz-zcell[i])*(midz-zcell[i])));

                axis.push_back(2);
                }

                //cout<<x<<" "<<y<<" "<<z<<" "<<t1<<" "<<t2<<" "<<xcell[i]<<" "<<ycell[i]<<" "<<zcell[i]<<" "<<sqrt((y-ycell[i])*(y-ycell[i])+(z-zcell[i])*(z-zcell[i]))<<" "<<i<<" codex"<<endl;

            }
    }


    if ((i > mod_straws[0]-1) && (i < mod_straws[1])) {
        if (y>= ycell[i] - 2*offset - outer && y<= ycell[i] + 2*offset + inner) {

            if (t1>0 && t2>0) {
                xhits.push_back(x);
                xcells.push_back(xcell[i]);

                yhits.push_back(y);
                ycells.push_back(ycell[i]);

                zhits.push_back(z);
                zcells.push_back(zcell[i]);

                radius.push_back(sqrt((midx-xcell[i])*(midx-xcell[i])+(midz-zcell[i])*(midz-zcell[i])));

                axis.push_back(1);
            }
                //cout<<x<<" "<<y<<" "    <<z<<" "<<t1<<" "<<t2<<" "<<xcell[i]<<" "<<ycell[i]<<" "<<zcell[i]<<" "<<sqrt((y-ycell[i])*(y-                                              ycell[i])+(z-zcell[i])*(z-zcell[i]))<<" "<<i<<" codex 2"<<endl;
            }
        }


    if ((i > mod_straws[2]-1) && (i < mod_straws[3])) {
        if (y>= ycell[i] - 2*offset - outer && y<= ycell[i] + 2*offset + inner) {

            if (t1>0 && t2>0) {
                xhits.push_back(x);
                xcells.push_back(xcell[i]);

                yhits.push_back(y);
                ycells.push_back(ycell[i]);

                zhits.push_back(z);
                zcells.push_back(zcell[i]);

                radius.push_back(sqrt((midx-xcell[i])*(midx-xcell[i])+(midz-zcell[i])*(midz-zcell[i])));

                axis.push_back(3);
            }
                //cout<<x<<" "<<y<<" "    <<z<<" "<<t1<<" "<<t2<<" "<<xcell[i]<<" "<<ycell[i]<<" "<<zcell[i]<<" "<<sqrt((y-ycell[i])*(y-                                              ycell[i])+(z-zcell[i])*(z-zcell[i]))<<" "<<i<<" codex 2"<<endl;
            }
    }

        }

    }

    if (xhits.size()<2) {

        xhits.push_back(0);
        xcells.push_back(0);

        yhits.push_back(0);
        ycells.push_back(0);

        zhits.push_back(0);
        zcells.push_back(0);
    }


    return std::make_tuple(xhits, yhits, zhits, xcells, ycells, zcells, radius, axis);

}

double linetocylinder(const double* params, Cylinder cylinder) {

    double x0 = params[0], y0 = params[1], z0 = params[2];  // Line origin
    double ux = params[3], uy = params[4], uz = params[5];  // Unit direction vector

    // Ensure the direction vector is normalized
    double norm = std::sqrt(ux * ux + uy * uy + uz * uz);
    ux /= norm; uy /= norm; uz /= norm;

    // Parametric line: p(t) = (x0 + t*ux, y0 + t*uy, z0 + t*uz)
    // Project the cylinder center onto the line
    double a = cylinder.a, b = cylinder.b, c = cylinder.c, r = cylinder.r, L = cylinder.L;
    int modu = cylinder.mod;

    double t = (a - x0) * ux + (b - y0) * uy + (c - z0) * uz;

    // Closest point on the line
    double px = x0 + t * ux;
    double py = y0 + t * uy;
    double pz = z0 + t * uz;

    double radial_distance;

    if (modu==0 || modu==2) {
        // Distance along the cylinder axis (x-axis)
        double dx = px - a;
        if (std::abs(dx) > L / 2) {
            // If out of cylinder bounds, penalize
            return 1e9;
        }

        // Perpendicular distance to the cylinder's surface
        double dy = py - b;
        double dz = pz - c;
        radial_distance = std::sqrt(dy * dy + dz * dz);

    }

    if (modu==1 || modu==3) {
        // Distance along the cylinder axis (x-axis)
        double dy = py - b;
        if (std::abs(dy) > L / 2) {
            // If out of cylinder bounds, penalize
            return 1e9;
        }

        // Perpendicular distance to the cylinder's surface
        double dx = px - b;
        double dz = pz - c;
        radial_distance = std::sqrt(dx * dx + dz * dz);

    }

    // Return the squared difference from the radius
    return (radial_distance - r) * (radial_distance - r);
}

double FitFunction(const double* params, const std::vector<Cylinder>& cylinders) {

    double x0 = params[0], y0 = params[1];
    double penalty = 0.0;

    // Add penalties for violating the constraints
    if (x0 < -0.3 || x0 > 0.3) penalty += 1e9 * (std::abs(x0) - 0.3);
    if (y0 < -0.3 || y0 > 0.3) penalty += 1e9 * (std::abs(y0) - 0.3);

    double sum = 0.0;
    for (const auto& cylinder : cylinders) {
        sum += linetocylinder(params, cylinder);
    }
    return sum + penalty;
}


std::pair<Point3D,Point3D> FitRadii(std::vector<Cylinder> cylinders, std::vector<double> first) {

    double initialParams[6] = {first[0], first[1], first[2], 0, 0, 1};
    auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    ROOT::Math::Functor func([&](const double* params) {
        return FitFunction(params, cylinders);}, 6);
    minimizer->SetFunction(func);

    // Set initial values and parameter step sizes
    minimizer->SetVariable(0, "x0", initialParams[0], 2);
    minimizer->SetVariable(1, "y0", initialParams[1], 2);
    minimizer->SetVariable(2, "z0", initialParams[2], 2);
    minimizer->SetVariable(3, "ux", initialParams[3], 0.1);
    minimizer->SetVariable(4, "uy", initialParams[4], 0.1);
    minimizer->SetVariable(5, "uz", initialParams[5], 0.1);

    // Perform the minimization
    minimizer->Minimize();

    // Retrieve results
    const double* results = minimizer->X();
    //std::cout << "Fitted line parameters:" << std::endl;

    //std::cout << "Origin: (" << results[0] << ", " << results[1] << ", " << results[2] << ")" << std::endl;
    //std::cout << "Direction: (" << results[3] << ", " << results[4] << ", " << results[5] << ")" << std::endl;

    Point3D centroid = {results[0],results[1],results[2]};
    Point3D direction = {results[3],results[4],results[5]};

    return std::make_pair(direction,centroid);

}

Point3D findClosestPointOnLine(const Point3D& centroid, const Point3D& direction, std::vector<double> point) {
    // Unpack components
    double xc = centroid[0], yc = centroid[1], zc = centroid[2];
    double dx = direction[0], dy = direction[1], dz = direction[2];
    double a = point[0], b = point[1], c = point[2];
    float dist = 99999;
    double tmin = 0;

    for (double t=-20; t<20; t=t+0.01) {
        float x = xc+t*dx;
        float y = yc+t*dy;
        float z = zc+t*dz;
        float diff = sqrt(pow(x-a,2)+pow(y-b,2)+pow(z-c,2));

        if (diff<dist) {
            tmin = t;
            dist=diff;
        }
    }

    // Compute the closest point
    Point3D closestPoint = {
        xc + tmin * dx, // x-coordinate
        yc + tmin * dy, // y-coordinate
        zc + tmin * dz  // z-coordinate
    };

    return closestPoint;
}

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);

  //std::ofstream outfile("shift_4mm_3layers_gap.csv", std::ios::app);
  std::ofstream outfile("dump", std::ios::app);

  //Generate straws
  auto [xCenter,yCenter,zCenter,mod_straws] = GetStrawCentersTransverseStepZoffsetOverhang();

  //Plot straw array
  //PlotXYCells(xCenter,yCenter,mod_straws);
  PlotMidYZ(zCenter,yCenter,mod_straws);

  //Generate and plot proton track
  //std::vector<double> trk_x,trk_y,trk_z,fit_x,fit_y,fit_z, dirvect=generateDir(), posvect=generatePos();
  //std::vector<double> trk_x,trk_y,trk_z,fit_x,fit_y,fit_z;//dirvect={0,1,-1}, posvect={0.1,0.1,0.1};

  int sample_N = 1000;
  double counts = 0;
  //for (double l = -7.5; l<7.5; l=l+0.1) {

      float sample = 0;
      float nhits = 0;
      float R_sum = 0;

      for (int k=0; k<sample_N; k++) {

    //std::vector<double> trk_x,trk_y,trk_z,dirvect = generateDir(), posvect = {0,0,0};

    auto [dirvect, posvect] = scatter();
    //posvect[2] = 0;

  //Determine hits
  auto [xhits, yhits, zhits, xcells, ycells, zcells, radius, axis] = hits(posvect, dirvect, xCenter, yCenter, zCenter, mod_straws);

  //cout<<xhits.size()<<endl;

  //Plot Hits, randomise DOCA 10%, plot DOCA
  //TGraph *graph = new TGraph();
  //graph->SetMarkerStyle(20); // Set marker style
  //graph->SetMarkerSize(1);   // Set marker size
  //graph->SetMarkerColor(3);   // Set marker size
  //graph->SetTitle("Example TGraph;X-axis;Y-axis");

  for (int i = 0; i < radius.size(); ++i) {
        radius[i] = randomiser(radius[i]);
        //graph->SetPoint(i, zhits[i], yhits[i]);
    }

    //graph->Draw("same, P");
    //PlotDOCAZY(zcells,ycells,radius);

    //Pack hits into cylinder objects and initalise track fitting algorithm

    std::vector<Cylinder> cylinders;
    std::vector<double> first = {0,0,1};
    double min = 9999;

    for (int i = 0; i<radius.size(); i++) {
        Cylinder cyl(xcells[i], ycells[i], zcells[i], radius[i], length, axis[i]);
        cylinders.push_back(cyl);

        double diff = sqrt(xcells[i]*xcells[i]+ycells[i]*ycells[i]+zcells[i]*zcells[i]);
        if (diff<min) {
            min = diff;
            first[0] = xcells[i];
            first[1] = ycells[i];
            first[2] = zcells[i];
        }

    }

    /*
    std::vector<double> trk_x,trk_y,trk_z;

    for (double t = 0; t<20; t=t+0.1) {
        trk_x.push_back(lab_dir[0]*t+position[0]);
        trk_y.push_back(lab_dir[1]*t+position[1]);
        trk_z.push_back(lab_dir[2]*t+position[2]);
    }

    PlotTrack(trk_z,trk_y);
    */


    //first = {0,0,1};

    //Fit cylinders
    auto [fitvec, fitcent] = FitRadii(cylinders,first);

    //Plot fitted track
    //for (double t = 0; t<20; t=t+0.1) {
        //fit_x.push_back(fitvec[0]*t+fitcent[0]);
        //fit_y.push_back(fitvec[1]*t+fitcent[1]);
        //fit_z.push_back(fitvec[2]*t+fitcent[2]);
    //}

    Point3D reco_v = findClosestPointOnLine(fitcent,fitvec,posvect);
    //PlotFit(fit_z,fit_y);

    if (xhits.size()>2) {
        //sample++;
        nhits=nhits+xhits.size();
        //reco_z->Fill(reco_v[2]-posvect[2]);
        //R_sum = R_sum + reco_v;
        reco_z->Fill(reco_v[2]-posvect[2]);

    }

    //reco_z_hits->Fill(xhits.size(),reco_v[2]-posvect[2]);

    if (xhits.size()>2) {
    counts = counts + 1;
    outfile<<"start"<<endl;
    outfile<<counts<<","<<posvect[0]<<","<<posvect[1]<<","<<posvect[2]<<","<<dirvect[0]<<","<<dirvect[1]<<","<<posvect[2]<<","<<reco_v[2]-posvect[2]<<","<<girth<<endl;

    for (int i=0; i<xhits.size(); i++) {
        outfile<<xhits[i]<<","<<yhits[i]<<","<<zhits[i]<<","<<radius[i]<<","<<axis[i]<<endl;
    }

    outfile<<"stop"<<endl;

    }



  }

     cout<<counts<<endl;
     cout<<counts/sample_N<<endl;
  //cout<<l<<","<<sample/3000.0<<","<<nhits/3000.0<<endl;

  //}

  //reco_z->Draw();

  //cout<<nhits/50000<<endl;;

  app.Run();
  return 0;
}
