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
#include <TPolyLine.h>

using namespace std;
using Point3D = array<double, 3>;
using Matrix3x3 = array<array<double, 3>, 3>;
const double c = 299792458.0;

float rCell = 0.5;
float vacuum = 150;
double offset = 35;
double length = 79.2;
int rows = 4;
double girth = 43.6;
float thick = 2;
float shift = 161;
float gap = 0.8;
double sig = 0;
int layers = 2;
int sample_N = 10000;

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

struct Strip {
    double a, b, c;  // Center position (a, b, c)
    double W;        // Radius
    double L;        // Length
    int mod;         // Module

    // Constructor
    Strip(double x, double y, double z, double width, double length, int stack)
        : a(x), b(y), c(z), W(width), L(length), mod(stack) {}
};

struct Point {
    double x, y, z;
};

double randomiser(double pos, double spread) {
    static std::random_device rd; // Random device for seeding
    static std::mt19937 gen(rd()); // Mersenne Twister generator
    static std::normal_distribution<> dis(0,spread);

    double randx = dis(gen);

    return randx;
}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>, std::vector<double>> GetBoxStrips() {

    std::vector<double> xCenter, yCenter, zCenter, mod;

    for (int k = 0; k<rows; k++) {

    for (int j = 0; j<layers; ++j) {
        xCenter.push_back(0);
        yCenter.push_back(length/2+j*thick+thick/2);
        zCenter.push_back(shift-(k+1)*girth-2*k*gap);
        mod.push_back(0);
    }

    for (int j = 0; j<layers; ++j) {
        xCenter.push_back(length/2+j*thick+thick/2);
        yCenter.push_back(0);
        zCenter.push_back(shift-(k+1)*girth-2*k*gap);
        mod.push_back(1);
    }

    for (int j = 0; j<layers; ++j) {
        xCenter.push_back(0);
        yCenter.push_back(-length/2-j*thick-thick/2);
        zCenter.push_back(shift-(k+1)*girth-2*k*gap);
        mod.push_back(2);
    }

    for (int j = 0; j<layers; ++j) {
        xCenter.push_back(-length/2-j*thick-thick/2);
        yCenter.push_back(0);
        zCenter.push_back(shift-(k+1)*girth-2*k*gap);
        mod.push_back(3);
    }

    }

    return std::make_tuple(xCenter, yCenter, zCenter, mod);

}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>, std::vector<double>> GetPentStrips() {

    std::vector<double> xCenter, yCenter, zCenter, mod;

    int sides = 5;

    for (int k = 0; k<rows; k++ ) {

        for (int j = 0; j<layers; ++j) {

            std::vector<float> vertx;
            std::vector<float> verty;

            for (int i = 0; i<sides; i++) {
                float circumradius = length/2+j*thick+thick/2;
                vertx.push_back(circumradius*cos(2*3.14159*i/sides));
                verty.push_back(circumradius*sin(2*3.14159*i/sides));
            }

            for (int i = 0; i<5; i++) {

                float midx = (vertx[i]+vertx[(i+1)])/2;
                float midy = (verty[i]+verty[(i+1)])/2;

                if (i==5) {
                    float midx = (vertx[i]+vertx[0])/2;
                    float midy = (verty[i]+verty[0])/2;
                }

                xCenter.push_back(midx);
                yCenter.push_back(midy);
                zCenter.push_back(shift-(k+1)*girth-2*k*gap);
                mod.push_back(1.0*i);
            }
        }
    }


    return std::make_tuple(xCenter, yCenter, zCenter, mod);

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

    TBox* cell = new TBox(-4.926, -3, 9.474, 3);
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

void PlotBoxYZ(const std::vector<double>& zCenter, const std::vector<double>& yCenter, std::vector<double> mod) {
    TCanvas* canvas = new TCanvas("canvas", "Straws and Track", 800, 800);
    // Create a canvas
    canvas->SetFixedAspectRatio(); // Ensure equal scaling on both axes
    canvas->DrawFrame(-200, -200, 200, 200); // Set the drawing frame (xMin, yMin, xMax, yMax)

    TEllipse* straw = new TEllipse(0,0,150);
            straw->SetFillStyle(0); // No fill
            straw->SetLineColor(kBlack);
            straw->Draw("same");

    TBox* cell = new TBox(-49.26, -30, 94.74, 30);
        cell->SetFillStyle(0); // No fill
        cell->SetLineColor(kRed);
        cell->Draw("same");

    TBox* vacuum = new TBox(-134, -63.5, 134, 63.5);
        vacuum->SetFillStyle(0); // No fill
        vacuum->SetLineColor(kBlue);
        vacuum->Draw("same");

    for (int i = 0; i < mod.size(); ++i) {

            if (mod[i]==0) {
                TBox* box1 = new TBox(zCenter[i]+girth/2, yCenter[i]+thick/2,zCenter[i]-girth/2,  yCenter[i]-thick/2);
                box1->SetFillStyle(0); // No fill
                box1->SetLineColor(kBlack);
                box1->Draw("same");
            }

            if (mod[i]==2) {
                TBox* box1 = new TBox(zCenter[i]+girth/2, yCenter[i]+thick/2,zCenter[i]-girth/2,  yCenter[i]-thick/2);
                box1->SetFillStyle(0); // No fill
                box1->SetLineColor(kBlack);
                box1->Draw("same");
            }


        }

    canvas->Update();
}

void PlotXYBox(const std::vector<double>& xCenter, const std::vector<double>& yCenter, std::vector<double> mod) {
    // Check if the inputs are valid

    TCanvas* canvas = new TCanvas("canvas", "Straws and Track", 1000, 1000);
    // Create a canvas
    canvas->SetFixedAspectRatio(); // Ensure equal scaling on both axes
    canvas->DrawFrame(-200, -200, 200, 200); // Set the drawing frame (xMin, yMin, xMax, yMax)

    TEllipse* DOCA = new TEllipse(0,0,vacuum);
            DOCA->SetFillStyle(0); // No fill
            DOCA->SetLineColor(kRed);
            DOCA->Draw("same");

    TEllipse* DOCA1 = new TEllipse(0,0,offset);
            DOCA1->SetFillStyle(0); // No fill
            DOCA1->SetLineColor(kBlack);
            DOCA1->Draw("same");


    // Draw the box (straws)

        for (int i = 0; i < mod.size(); ++i) {

            if (mod[i]==0) {
                TBox* box1 = new TBox(xCenter[i]+length/2, yCenter[i]+thick/2,xCenter[i]-length/2,  yCenter[i]-thick/2);
                box1->SetFillStyle(0); // No fill
                box1->SetLineColor(kBlack);
                box1->Draw("same");
            }

            if (mod[i]==1) {
                TBox* box1 = new TBox(xCenter[i]+thick/2, yCenter[i]+length/2,xCenter[i]-thick/2,  yCenter[i]-length/2);
                box1->SetFillStyle(0); // No fill
                box1->SetLineColor(kBlack);
                box1->Draw("same");
            }

            if (mod[i]==2) {
                TBox* box1 = new TBox(xCenter[i]+length/2, yCenter[i]+thick/2,xCenter[i]-length/2,  yCenter[i]-thick/2);
                box1->SetFillStyle(0); // No fill
                box1->SetLineColor(kBlack);
                box1->Draw("same");
            }

            if (mod[i]==3) {
                TBox* box1 = new TBox(xCenter[i]-thick/2, yCenter[i]+length/2,xCenter[i]+thick/2,  yCenter[i]-length/2);
                box1->SetFillStyle(0); // No fill
                box1->SetLineColor(kBlack);
                box1->Draw("same");
            }


        }
    /*

    for (int i = mod_straws[0]; i < mod_straws[1]; ++i) {
        TBox* straw1 = new TBox(xCenter[i]-rCell, yCenter[i]-girth/2, xCenter[i]+rCell, yCenter[i]+girth/2-3);
        straw1->SetFillStyle(0); // No fill
        straw1->SetLineColor(kRed);
        straw1->Draw("same");
    }

    for (int i = mod_straws[1]; i < mod_straws[2]; ++i) {
        TBox* straw2 = new TBox(xCenter[i]-girth/2, yCenter[i]-rCell, xCenter[i]+girth/2, yCenter[i]+rCell);
        straw2->SetFillStyle(0); // No fill
        straw2->SetLineColor(kBlue);
        straw2->Draw("same");
    }

    for (int i = mod_straws[2]; i < mod_straws[3]; ++i) {
        TBox* straw3 = new TBox(xCenter[i]-rCell, yCenter[i]-girth/2, xCenter[i]+rCell, yCenter[i]+girth/2);
        straw3->SetFillStyle(0); // No fill
        straw3->SetLineColor(kGreen);
        straw3->Draw("same");

    }

    */

    // Update the canvas to render the drawing
    canvas->Update();
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>,
           std::vector<double>, std::vector<double>, std::vector<int>>
BoxHits(std::vector<double> posvect, std::vector<double> dirvect, std::vector<double> xcell, std::vector<double> ycell, std::vector<double> zcell, std::vector<double> mod) {

    std::vector<double> xcells, ycells, zcells;
    std::vector<double> xhits, yhits, zhits;
    std::vector<int> axis;

    double mag = sqrt(dirvect[0]*dirvect[0]+dirvect[1]*dirvect[1]+dirvect[2]*dirvect[2]);
    double ux = dirvect[0]/mag, uy = dirvect[1]/mag, uz = dirvect[2]/mag;
    double x0 = posvect[0], y0 = posvect[1], z0 = posvect[2];
    float L = length, W = girth, T = thick;

    for (int i=0; i<xcell.size(); i++) {

    double rx = xcell[i], ry = ycell[i], rz = zcell[i];

    L = thick/2, W = girth, T = length;

    if (mod[i] == 0) {

    // Intersection with x = rx + T/2
    if (ux != 0) {
        double t = (rx + T/2 - x0) / ux;
        double y = y0 + t * uy;
        double z = z0 + t * uz;
        if (ry - L/2 <= y && y <= ry + L/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(rx+T/2);
            yhits.push_back(y);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(0);
            continue;
        }
    }

    // Intersection with x = rx - T/2
    if (ux != 0) {
        double t = (rx + T/2 - x0) / ux;
        double y = y0 + t * uy;
        double z = z0 + t * uz;
        if (ry - L/2 <= y && y <= ry + L/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(rx-T/2);
            yhits.push_back(y);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(0);
            continue;
        }
    }

    // Intersection with y = ry - L/2
    if (uy != 0) {
        double t = (ry - L/2 - y0) / uy;
        double x = x0 + t * ux;
        double z = z0 + t * uz;
        if (rx - T/2 <= x && x <= rx + T/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(ry-L/2);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(0);
            continue;
        }
    }

    // Intersection with y = ry + L/2
    if (uy != 0) {
        double t = (ry + L/2 - y0) / uy;
        double x = x0 + t * ux;
        double z = z0 + t * uz;
        if (rx - T/2 <= x && x <= rx + T/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(ry+L/2);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(0);
            continue;
        }
    }

    // Intersection with z = rz - W/2
    if (uz != 0) {
        double t = (rz - W/2 - z0) / uz;
        double x = x0 + t * ux;
        double y = y0 + t * uy;
        if (rx - T/2 <= x && x <= rx + T/2 && ry - L/2 <= y && y <= ry + L/2 && t>0) {
           xhits.push_back(x);
            yhits.push_back(y);
            zhits.push_back(rz-W/2);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(0);
            continue;
        }
    }

    // Intersection with z = rz + W/2
    if (uz != 0) {
        double t = (rz + W/2 - z0) / uz;
        double x = x0 + t * ux;
        double y = y0 + t * uy;
        if (rx - T/2 <= x && x <= rx + T/2 && ry - L/2 <= y && y <= ry + L/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(y);
            zhits.push_back(rz+W/2);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(0);
            continue;
        }
    }

    }

    L = length, W = girth, T = thick/2;

    if (mod[i] == 1) {

    // Intersection with x = rx - T/2
    if (ux != 0) {
        double t = (rx - T/2 - x0) / ux;
        double y = y0 + t * uy;
        double z = z0 + t * uz;
        if (ry - L/2 <= y && y <= ry + L/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(rx-T/2);
            yhits.push_back(y);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(1);
            continue;
        }
    }

    // Intersection with x = rx + T/2
    if (ux != 0) {
        double t = (rx + T/2 - x0) / ux;
        double y = y0 + t * uy;
        double z = z0 + t * uz;
        if (ry - L/2 <= y && y <= ry + L/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(rx+T/2);
            yhits.push_back(y);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(1);
            continue;
        }
    }

    // Intersection with y = ry - L/2
    if (uy != 0) {
        double t = (ry - L/2 - y0) / uy;
        double x = x0 + t * ux;
        double z = z0 + t * uz;
        if (rx - T/2 <= x && x <= rx + T/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(ry-L/2);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(1);
            continue;
        }
    }

    // Intersection with y = ry + L/2
    if (uy != 0) {
        double t = (ry + L/2 - y0) / uy;
        double x = x0 + t * ux;
        double z = z0 + t * uz;
        if (rx - T/2 <= x && x <= rx + T/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(ry+L/2);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(1);
            continue;
        }
    }

    // Intersection with z = rz - W/2
    if (uz != 0) {
        double t = (rz - W/2 - z0) / uz;
        double x = x0 + t * ux;
        double y = y0 + t * uy;
        if (rx - T/2 <= x && x <= rx + T/2 && ry - L/2 <= y && y <= ry + L/2 && t>0) {
           xhits.push_back(x);
            yhits.push_back(y);
            zhits.push_back(rz-W/2);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(1);
            continue;
        }
    }

    // Intersection with z = rz + W/2
    if (uz != 0) {
        double t = (rz + W/2 - z0) / uz;
        double x = x0 + t * ux;
        double y = y0 + t * uy;
        if (rx - T/2 <= x && x <= rx + T/2 && ry - L/2 <= y && y <= ry + L/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(y);
            zhits.push_back(rz+W/2);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(1);
            continue;
        }
    }

    }

    L = thick/2, W = girth, T = length;

    if (mod[i] == 2) {

    // Intersection with x = rx - T/2
    if (ux != 0) {
        double t = (rx - T/2 - x0) / ux;
        double y = y0 + t * uy;
        double z = z0 + t * uz;
        if (ry - L/2 <= y && y <= ry + L/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(rx-T/2);
            yhits.push_back(y);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(2);
            continue;
        }
    }

    // Intersection with x = rx + T/2
    if (ux != 0) {
        double t = (rx + T/2 - x0) / ux;
        double y = y0 + t * uy;
        double z = z0 + t * uz;
        if (ry - L/2 <= y && y <= ry + L/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(rx+T/2);
            yhits.push_back(y);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(2);
            continue;
        }
    }

    // Intersection with y = ry - L/2
    if (uy != 0) {
        double t = (ry - L/2 - y0) / uy;
        double x = x0 + t * ux;
        double z = z0 + t * uz;
        if (rx - T/2 <= x && x <= rx + T/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(ry-L/2);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(2);
            continue;
        }
    }

    // Intersection with y = ry + L/2
    if (uy != 0) {
        double t = (ry + L/2 - y0) / uy;
        double x = x0 + t * ux;
        double z = z0 + t * uz;
        if (rx - T/2 <= x && x <= rx + T/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(ry+L/2);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(2);
            continue;
        }
    }

    // Intersection with z = rz - W/2
    if (uz != 0) {
        double t = (rz - W/2 - z0) / uz;
        double x = x0 + t * ux;
        double y = y0 + t * uy;
        if (rx - T/2 <= x && x <= rx + T/2 && ry - L/2 <= y && y <= ry + L/2 && t>0) {
           xhits.push_back(x);
            yhits.push_back(y);
            zhits.push_back(rz-W/2);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(2);
            continue;
        }
    }

    // Intersection with z = rz + W/2
    if (uz != 0) {
        double t = (rz + W/2 - z0) / uz;
        double x = x0 + t * ux;
        double y = y0 + t * uy;
        if (rx - T/2 <= x && x <= rx + T/2 && ry - L/2 <= y && y <= ry + L/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(y);
            zhits.push_back(rz+W/2);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(2);
            continue;
        }
    }

    }

    L = length, W = girth, T = thick/2;

    if (mod[i] == 3) {

    // Intersection with x = rx - T/2
    if (ux != 0) {
        double t = (rx - T/2 - x0) / ux;
        double y = y0 + t * uy;
        double z = z0 + t * uz;
        if (ry - L/2 <= y && y <= ry + L/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(rx-T/2);
            yhits.push_back(y);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(3);
            continue;
        }
    }

    // Intersection with x = rx + T/2
    if (ux != 0) {
        double t = (rx + T/2 - x0) / ux;
        double y = y0 + t * uy;
        double z = z0 + t * uz;
        if (ry - L/2 <= y && y <= ry + L/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(rx+T/2);
            yhits.push_back(y);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(3);
            continue;
        }
    }

    // Intersection with y = ry - L/2
    if (uy != 0) {
        double t = (ry - L/2 - y0) / uy;
        double x = x0 + t * ux;
        double z = z0 + t * uz;
        if (rx - T/2 <= x && x <= rx + T/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(ry-L/2);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(3);
            continue;
        }
    }

    // Intersection with y = ry + L/2
    if (uy != 0) {
        double t = (ry + L/2 - y0) / uy;
        double x = x0 + t * ux;
        double z = z0 + t * uz;
        if (rx - T/2 <= x && x <= rx + T/2 && rz - W/2 <= z && z <= rz + W/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(ry+L/2);
            zhits.push_back(z);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(3);
            continue;
        }
    }

    // Intersection with z = rz - W/2
    if (uz != 0) {
        double t = (rz - W/2 - z0) / uz;
        double x = x0 + t * ux;
        double y = y0 + t * uy;
        if (rx - T/2 <= x && x <= rx + T/2 && ry - L/2 <= y && y <= ry + L/2 && t>0) {
           xhits.push_back(x);
            yhits.push_back(y);
            zhits.push_back(rz-W/2);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(3);
            continue;
        }
    }

    // Intersection with z = rz + W/2
    if (uz != 0) {
        double t = (rz + W/2 - z0) / uz;
        double x = x0 + t * ux;
        double y = y0 + t * uy;
        if (rx - T/2 <= x && x <= rx + T/2 && ry - L/2 <= y && y <= ry + L/2 && t>0) {
            xhits.push_back(x);
            yhits.push_back(y);
            zhits.push_back(rz+W/2);
            xcells.push_back(xcell[i]);
            ycells.push_back(ycell[i]);
            zcells.push_back(zcell[i]);
            axis.push_back(3);
            continue;
        }
    }

    }

    L = thick/2, W = girth, T = length;

    }

    if (xhits.size()<2) {

        xhits.push_back(0);
        xcells.push_back(0);

        yhits.push_back(0);
        ycells.push_back(0);

        zhits.push_back(0);
        zcells.push_back(0);
    }


    return std::make_tuple(xhits, yhits, zhits, xcells, ycells, zcells, axis);

}

void PlotXYPent(const std::vector<double>& xCenter, const std::vector<double>& yCenter, std::vector<double> mod) {
    TCanvas* canvas = new TCanvas("canvas", "Straws and Track", 1000, 1000);
    canvas->SetFixedAspectRatio();
    canvas->DrawFrame(-200, -200, 200, 200);

    // Draw reference circles
    TEllipse* DOCA = new TEllipse(0, 0, vacuum);
    DOCA->SetFillStyle(0);
    DOCA->SetLineColor(kRed);
    DOCA->Draw("same");

    TEllipse* DOCA1 = new TEllipse(0, 0, offset);
    DOCA1->SetFillStyle(0);
    DOCA1->SetLineColor(kBlack);
    DOCA1->Draw("same");

    const int sides = 5; // Pentagon

    for (int j = 0; j < layers; j++) {
        double vertx[sides + 1]; // Ensure size is 6
        double verty[sides + 1];

        double circumradius = length/(2*sin(M_PI/5))+j*thick;

        for (int i = 0; i < sides; i++) {
            double angle = 2 * M_PI * i / sides;
            vertx[i] = circumradius * cos(angle);
            verty[i] = circumradius * sin(angle);
        }

        // Close the pentagon by repeating the first vertex
        vertx[sides] = vertx[0];
        verty[sides] = verty[0];

        // Draw the pentagon using TPolyLine
        TPolyLine* pentagon = new TPolyLine(sides + 1, vertx, verty);
        pentagon->SetLineColor(kBlue);
        pentagon->SetFillStyle(0);
        pentagon->SetLineWidth(1);
        pentagon->Draw("same");
    }

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
  graph->SetMarkerStyle(3);  // 21: circle
  graph->SetMarkerColor(kBlue);  // Red points
  graph->SetLineStyle(0);  // Red points
  graph->SetMarkerSize(0.3);  // Marker size

  // Draw the graph on the canvas
  graph->Draw("same, P");

}

void PlotNoise(std::vector<double> trk_x,std::vector<double> trk_y) {

  TGraph* graph = new TGraph(trk_x.size(), trk_x.data(), trk_y.data());

  // Set the marker style for the points (e.g., a circle)
  graph->SetMarkerStyle(3);  // 21: circle
  graph->SetMarkerColor(kRed);  // Red points
  graph->SetLineStyle(0);  // Red points
  graph->SetMarkerSize(0.3);  // Marker size

  // Draw the graph on the canvas
  graph->Draw("same, P");

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
    static std::normal_distribution<> disXY(0,2);//XY dispursion
    static std::uniform_real_distribution<> disAB(-5.0/1000.0,5.0/1000.0);//Angular dispursion
    //static std::uniform_real_distribution<> disZ(-49.26,0.74);//Small target z
    static std::uniform_real_distribution<> disZ(-49.26,94.74);//Long target z


    // Generate incoming angle and position
    double x0 = disXY(gen), y0 = disXY(gen);
    double A = disAB(gen), B = disAB(gen);
    double z0 = sqrt(1.0 - tan(A)*tan(A) - tan(B)*tan(B));
    std::vector<double> direction = {tan(A),tan(B),z0};

    //Calculate incoming direction and positions
    //double z = disZ(gen) + 57.5;//small target
    double z = disZ(gen);//large target
    double x = x0 + z*tan(A), y = y0 + z*tan(B);
    std::vector<double> position = {x,y,z};

    return position;

}

double lorentzFactor(double v) {
    return 1.0 / std::sqrt(1.0 - (v * v));
}

std::vector<double> lorentzBoostDirection(const std::vector<double>& particleDirection, double sourceVelocity) {
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

     std::vector<double> trk_x,trk_y,trk_z;

    for (double t = 0; t<20; t=t+0.1) {
        trk_x.push_back(lab_dir[0]*t+position[0]);
        trk_y.push_back(lab_dir[1]*t+position[1]);
        trk_z.push_back(lab_dir[2]*t+position[2]);
    }

    //PlotTrack(trk_z,trk_y);

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

std::pair<float,float> randcirc() {
    // Seed the random number generator
    static bool seeded = false;
    if (!seeded) {
        std::srand(static_cast<unsigned int>(std::time(0)));
        seeded = true;
    }

    // Generate a random angle between 0 and 2π
    double angle = 2.0 * M_PI * static_cast<double>(std::rand()) / RAND_MAX;

    // Convert the angle to Cartesian coordinates

    return std::make_pair(cos(angle),sin(angle));
}


int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);

  std::ofstream outfile("box4.csv", std::ios::app);

  //Generate straws
  auto [xCenter,yCenter,zCenter,mod] = GetBoxStrips();
  //auto [xCenter,yCenter,zCenter,mod] = GetPentStrips();
  int counts = 0;

  //Plot straw array
  //PlotXYBox(xCenter,yCenter,mod);
  //PlotXYPent(xCenter,yCenter,mod);
  PlotBoxYZ(zCenter,yCenter,mod);

  for (int l = 0; l<sample_N; l++) {

  //Generate and plot proton track
  std::vector<double> trk_x,trk_y,trk_z; //dirvect={0,0,0}, posvect={0.0,0.0,0.0};

  //auto [x,y] = randcirc();

  //dirvect={x,y,0};

   auto [dirvect, posvect] = scatter();
  //posvect={0.1,0.1,0.1};
  //rvect[0] = 0;
  //dirvect={0,1,0};

  for (double i= 0; i<200.0; i = i + 0.5) {
      trk_x.push_back(posvect[0]+i*dirvect[0]);
      trk_y.push_back(posvect[1]+i*dirvect[1]);
      trk_z.push_back(posvect[2]+i*dirvect[2]);
  }

  //Determine hits
  auto [xhits, yhits, zhits, xcells, ycells, zcells, axis] = BoxHits(posvect, dirvect, xCenter, yCenter, zCenter, mod);

  PlotTrack(trk_z, trk_y);
  PlotFit(zhits,yhits);

  for (int i=0; i<xhits.size(); i++) {
        cout<<xhits[i]<<","<<yhits[i]<<","<<zhits[i]<<","<<xcells[i]<<","<<ycells[i]<<","<<zcells[i]<<endl;
  }

  for (int i = 0; i < xhits.size(); ++i) {

        float rand = randomiser(zhits[i],1.2/2.35);

    }

  //PlotNoise(xhits,yhits);

    if (xhits.size()>1) {
    counts = counts + 1;
    outfile<<"start"<<endl;
    outfile<<counts<<","<<sample_N<<","<<posvect[0]<<","<<posvect[1]<<","<<posvect[2]<<","<<dirvect[0]<<","<<dirvect[1]<<","<<posvect[2]<<endl;

    for (int i=0; i<xhits.size(); i++) {
        outfile<<xhits[i]<<","<<yhits[i]<<","<<zhits[i]<<","<<length<<","<<girth<<","<<axis[i]<<endl;
    }

    outfile<<"stop"<<endl;

    }

     cout<<counts<<endl;
     cout<<counts/sample_N<<endl;


  //cout<<l<<","<<sample/3000.0<<","<<nhits/3000.0<<endl;

  }

  //reco_z->Draw();

  //cout<<nhits/50000<<endl;;

  app.Run();
  return 0;
}
