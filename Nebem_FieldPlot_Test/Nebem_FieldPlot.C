#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>
#include "TString.h"

// ROOT includes
#include <TFile.h>
#include <TH2D.h>
#include <TGeoManager.h>
#include <TTree.h>

#include <ctime>
#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
// Standard includes

#include <iomanip>
#include <stdexcept>
#include <fstream>
#include <array>

#include "Garfield/AvalancheGrid.hh"
#include "Garfield/AvalancheMicroscopic.hh"
//#include "Garfield/ComponentParallelPlate.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewGeometry.hh"

// Setup geometry and material
#include "Garfield/ComponentNeBem3d.hh"
#include "Garfield/GeometrySimple.hh"
#include <Garfield/SolidBox.hh>
#include <Garfield/Solid.hh>

// Medium Classes
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Medium.hh"
#include "Garfield/MediumConductor.hh"
#include "Garfield/MediumPlastic.hh"
#include "Garfield/FundamentalConstants.hh"

#define LOG(x) std::cout << x << std::endl

using namespace Garfield;

int main(int argc, char *argv[]) {
  // Creates ROOT Application Environtment
  TApplication app("app", &argc, argv);

  const bool debug = true;
  constexpr bool plotSignal = true;
  constexpr bool plotField = true;
  
   static constexpr double fGasGapThickness = 0.025;         // cm - total gas gap thickness            250 micron = 0.25 mm = 0.025 cm
    static constexpr double fGasGapCenterY = 0.0/10.0;           // cm - centered at Y=0                    y = 0 mm =0 cm
    static constexpr double fAnodeCathodeThickness = 0.02/10.0;  // cm - thickness of HV layers             20 micron = 0.02 mm = 0.002 cm
    static constexpr double fStripThickness = 0.05/10.0;         // cm - thickness of copper strips         50 micron = 0.005 cm
    static constexpr double fStripWidthX = 30.0/10;              // cm - readout strip width                30 mm = 3 cm
    static constexpr double fDetectorSizeX = 30.0/10.0;         // cm - detector size in X                 30 mm = 3 cm    
    static constexpr double fDetectorSizeZ = 30.0/10.0;         // cm - detector size in Z                 30 mm = 3 cm
    static constexpr double fHoneyCombThickness = 0;     // cm - honey comb layer thickness         0 mm = 0 cm (Remove Honeycomb)
    static constexpr double fMylarThickness = 0.035;         // cm - mylar layer thickness              350 micron = 0.35 mm = 0.035 cm
    static constexpr double fResistiveGlassThickness = 0.12; // cm - resistive glass layer thickness  1.2 mm = 0.12 cm
    static constexpr double fAnodeVoltage = 5500.0;               // V - ANODE at +5.5 kV
    static constexpr double fCathodeVoltage = -5500.0;            // V - CATHODE at -5.5 kV
    
    // Calculate the dimension and location of geometry.
    double gasGapTop = fGasGapCenterY + fGasGapThickness/2.0;
    double gasGapBottom = fGasGapCenterY - fGasGapThickness/2.0;
    double fResistiveGlassPosition = fGasGapCenterY + fGasGapThickness/2.0 + fResistiveGlassThickness/2.0;
    double fAnodeCathodePosition = fGasGapCenterY + fGasGapThickness/2.0 + fResistiveGlassThickness + fAnodeCathodeThickness/2.0;
    double fMylarPosition = fGasGapCenterY + fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness/2.0;
    double fStripPosition = fGasGapCenterY + fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness + fStripThickness/2.0;

  
MediumMagboltz gas;
gas.SetComposition("C2H2F4",90.0,"iC4H10",5,"SF6",5);
gas.SetTemperature(296.15); 
gas.SetPressure(760.0);
gas.Initialise(true);

  // Materials needed.
    MediumPlastic* glass = new Garfield::MediumPlastic();
    MediumConductor* graphite = new Garfield::MediumConductor();
    MediumPlastic* mylar = new Garfield::MediumPlastic();
    MediumConductor* copper = new Garfield::MediumConductor();

    glass->SetDielectricConstant(10.0);
    mylar->SetDielectricConstant(3.1);

    // Define Geometry
    GeometrySimple* rpc_geometry = new Garfield::GeometrySimple();

    // Define Component
    ComponentNeBem3d *rpc = new ComponentNeBem3d();

    // Readout Strips
    Garfield::SolidBox* Box_TopStrip = new Garfield::SolidBox(0,  fStripPosition ,
                                                               0, fDetectorSizeX/2,
                                                               fStripThickness/2, fDetectorSizeZ/2);
    Box_TopStrip->SetLabel("TopReadout");
    Box_TopStrip->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_TopStrip,copper);
    const std::string label = "TopReadout";

    Garfield::SolidBox* Box_BottomStrip = new Garfield::SolidBox(0, - fStripPosition,
                                                               0, fDetectorSizeX/2,
                                                               fStripThickness/2, fDetectorSizeZ/2);
    Box_BottomStrip->SetLabel("BottomReadout");
    Box_BottomStrip->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_BottomStrip,copper);
    
    // Mylar Layers
    Garfield::SolidBox* Box_TopMylar = new Garfield::SolidBox(0, fMylarPosition,
                                                               0, fDetectorSizeX/2,
                                                               fMylarThickness/2, fDetectorSizeZ/2);
    Box_TopMylar->SetLabel("TopMylar");
    Box_TopMylar->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_TopMylar,mylar);

    Garfield::SolidBox* Box_BottomMylar = new Garfield::SolidBox(0, -fMylarPosition,
                                                               0, fDetectorSizeX/2,
                                                               fMylarThickness/2, fDetectorSizeZ/2);
    Box_BottomMylar->SetLabel("BottomMylar");
    Box_BottomMylar->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_BottomMylar,mylar);
    
    // Graphite Layers (Anode and Cathode)
    Garfield::SolidBox* Box_TopGraphite = new Garfield::SolidBox(0, fAnodeCathodePosition,
                                                               0, fDetectorSizeX/2,
                                                               fAnodeCathodeThickness/2, fDetectorSizeZ/2);
    Box_TopGraphite->SetLabel("TopGraphite");
    Box_TopGraphite->SetBoundaryPotential(fAnodeVoltage); // Anode
    rpc_geometry->AddSolid(Box_TopGraphite,graphite);

    Garfield::SolidBox* Box_BottomGraphite = new Garfield::SolidBox(0, -fAnodeCathodePosition,
                                                               0, fDetectorSizeX/2,
                                                               fAnodeCathodeThickness/2, fDetectorSizeZ/2);
    Box_BottomGraphite->SetLabel("BottomGraphite");
    Box_BottomGraphite->SetBoundaryPotential(fCathodeVoltage); // Cathode
    rpc_geometry->AddSolid(Box_BottomGraphite,graphite);

    // Resistive Glass Layers
    Garfield::SolidBox* Box_TopResistiveGlass = new Garfield::SolidBox(0, fResistiveGlassPosition,
                                                               0, fDetectorSizeX/2,
                                                               fResistiveGlassThickness/2, fDetectorSizeZ/2);
    Box_TopResistiveGlass->SetLabel("TopResistiveGlass");
    Box_TopResistiveGlass->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_TopResistiveGlass,glass);
    Garfield::SolidBox* Box_BottomResistiveGlass = new Garfield::SolidBox(0, -fResistiveGlassPosition,
                                                               0, fDetectorSizeX/2,
                                                               fResistiveGlassThickness/2, fDetectorSizeZ/2);
    Box_BottomResistiveGlass->SetLabel("BottomResistiveGlass");
    Box_BottomResistiveGlass->SetBoundaryDielectric();
    rpc_geometry->AddSolid(Box_BottomResistiveGlass,glass);

    // Gas Gap
    // Add gas gap volume
    Garfield::SolidBox* Box_GasGap = new Garfield::SolidBox(
        0, fGasGapCenterY, 0,                  // center
        fDetectorSizeX/2,                      // half-width X
        fGasGapThickness/2,                    // half-width Y
        fDetectorSizeZ/2);
    Box_GasGap->SetLabel("GasGap");
    rpc_geometry->AddSolid(Box_GasGap, &gas);           // Use the gas medium

    rpc->SetGeometry(rpc_geometry);
    rpc->SetTargetElementSize(0.01); // 0.01 cm = 0.1mm = 100 microns  
    rpc->SetMinMaxNumberOfElements(3, 10);  //
    rpc->EnableDebugging();
    rpc->Initialise();

  // Create the sensor.
  Sensor sensor(rpc);
  sensor.AddElectrode(rpc, label);
  
TCanvas* cGeo = new TCanvas("c","Geometry",800,600);
ViewGeometry geoView(rpc_geometry);
geoView.SetCanvas(cGeo);
geoView.Plot3d();
cGeo->Update();
cGeo->SaveAs("../GeometryPlot.png");

const std::size_t nx = 50,ny =50,ncont =104;
const double dTotal = 2*(fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness + fStripThickness);


  ViewField *fieldView = nullptr;
  TCanvas *cField = nullptr;
  cField = new TCanvas("cField", "Field", 1200,800);
  fieldView = new ViewField();
  fieldView->SetCanvas(cField);
  fieldView->SetComponent(rpc);
  
  fieldView->SetNumberOfSamples2d(nx,ny);
  fieldView->Plot("emag");
  
 //fieldView->EnableAutoRange();
 fieldView->SetElectricFieldRange(0,10000); 
 fieldView->PlotProfile(0,-0.025,0,0,0.025,0,"emag",true);


  ViewField *contourView = nullptr;
  TCanvas *cContour = nullptr;
  cContour = new TCanvas("cContour","Contour",1200,800);
  contourView = new ViewField();
  contourView->SetCanvas(cContour);
  contourView->SetComponent(rpc);
  contourView->SetPlane(0,-1,0,0,0.0125,0);
  contourView->SetArea(-5,-0.025,-5,5,0.02,5);
  contourView->SetNumberOfContours(ncont);
  contourView->PlotContour("emag");
  cContour->SaveAs("../contouronxz.png");

  ViewField *field2View = nullptr;
  TCanvas *cfield2 = nullptr;
  cfield2 = new TCanvas("cfield2","field2",1200,800);
  field2View = new ViewField();
  field2View->SetCanvas(cfield2);
  field2View->SetComponent(rpc);
  field2View->SetPlane(0,-1,0,0,0.015,0);
  field2View->SetArea(-5,-0.025,-5,5,0.02,5);
  field2View->SetNumberOfSamples2d(nx,ny);
  field2View->Plot("emag");
  cfield2->SaveAs("../fieldonxz.png");

struct fielddata
{
  int pointindex;
  double x;
  double y;
  double z;
  double ex;
  double ey;
  double ez;
  double emag;
  int status;
};

fielddata NebemField;

// Construct proper filename
std::string filename = __FILE_NAME__;
auto pos = filename.find_last_of('.');
if (pos != std::string::npos) {
    filename = filename.substr(0, pos);
}
filename += ".root";

TFile *treefile = new TFile(filename.c_str(),"recreate");

TTree *fieldtree = new TTree("FieldTree","Field value in gas gap.");

fieldtree->Branch("point index",&NebemField.pointindex,"pointindex/I");
fieldtree->Branch("x",&NebemField.x,"x/D");
fieldtree->Branch("y",&NebemField.y,"y/D");
fieldtree->Branch("z",&NebemField.z,"z/D");
fieldtree->Branch("ex",&NebemField.ex,"ex/D");
fieldtree->Branch("ey",&NebemField.ey,"ey/D");
fieldtree->Branch("ez",&NebemField.ez,"ez/D");
fieldtree->Branch("emag",&NebemField.emag,"emag/D");
fieldtree->Branch("status",&NebemField.status,"status/I");

double xmin = -1.5, xmax = 1.5;
double zmin = -1.5, zmax = 1.5;

double ymin = -fGasGapThickness/2 , ymax = +fGasGapThickness;

int nfx = 50, nfz = 50;
int nfy = 50;

double dx = (xmax - xmin)/(double)(nfx-1);
double dz = (zmax - zmin)/(double)(nfz-1);
double dy = (ymax - ymin)/(double)(nfy-1);

int index = 0;
for (int i=0; i<nfx; ++i) // loop over all x grid points
{
  double x = xmin + i*dx;
  for (int j=0; j<nfy; ++j)
  {
    double y = ymin + j*dy;
    for (int k=0; k<nfz; ++k)
    {
      double z = zmin + k*dz;
      Medium* medium = nullptr;
      rpc->ElectricField(x,y,z,NebemField.ex,NebemField.ey,NebemField.ez,medium,NebemField.status);  
      if (NebemField.status == 0) // field in gas medium.
      {
        NebemField.x = x;
        NebemField.y = y;
        NebemField.z = z;
        NebemField.emag = std::sqrt(NebemField.ex*NebemField.ex + NebemField.ey*NebemField.ey + NebemField.ez*NebemField.ez);
        NebemField.pointindex = index++;
        fieldtree->Fill();

        if (index%100 ==0){
          std::cout<<"Sampled " << index << " points." << std::endl;
        }
      }
    }
  }
}
std::cout << gDirectory->GetName() << std::endl;
fieldtree->Write();
treefile->Close();

  /*
  // Set the time bins.
  const std::size_t nTimeBins = 200;
  const double tmin = 0.;
  const double tmax = 4;
  const double tstep = (tmax - tmin) / nTimeBins;
  sensor.SetTimeWindow(tmin, tstep, nTimeBins);

  // Use microscopic tracking for initial stage of the avalanche.
  AvalancheMicroscopic aval(&sensor);
  // Set time until which the calculations will be done microscopically.
  const double tMaxWindow = 4;
  aval.SetTimeWindow(0., tMaxWindow);

  // Use a grid-based method for simulating the avalanche growth 
  // after the initial stage.
  //AvalancheGrid avalgrid(&sensor);
  
  // Start the track in the first gas layer.
  const double dTotal = 2*(fGasGapThickness/2.0 + fAnodeCathodeThickness + fResistiveGlassThickness + fMylarThickness + fStripThickness);
  const int nX = 5;
  const int nZ = 5;
  const double dY = 1.e-4;
  const int nY = int(dTotal / dY);
  //avalgrid.SetGrid(-0.05, 0.05, nX, -dTotal/2, dTotal/2, nY, -0.05, 0.05, nZ);

  // Prepare the plotting of the induced current and charge.
  ViewSignal *signalView = nullptr;
  TCanvas *cSignal = nullptr;
  if (plotSignal) {
    cSignal = new TCanvas("cSignal", "Current - 1 event", 600, 600);
    signalView = new ViewSignal(&sensor);
    signalView->SetCanvas(cSignal);
  }

  ViewSignal *chargeView = nullptr;
  TCanvas *cCharge = nullptr;
  if (plotSignal) {
    cCharge = new TCanvas("cCharge", "Charge - 1 event", 600, 600);
    chargeView = new ViewSignal(&sensor);
    chargeView->SetCanvas(cCharge);
  }

  // Set up Heed.
  TrackHeed track(&sensor);
  // Set the particle type and momentum [eV/c].
  track.SetParticle("pion");
  track.SetMomentum(7.e9);
  track.CrossInactiveMedia(true);

  std::clock_t start = std::clock();

  auto cD = new TCanvas("cD","Drift",600,600);
  ViewDrift driftView;
  driftView.SetCanvas(cD);
  aval.EnablePlotting(&driftView, 100);
  track.EnablePlotting(&driftView);


  const double y0 = dTotal/2 - fStripThickness - fMylarThickness - fAnodeCathodeThickness - fResistiveGlassThickness;
  // Simulate a charged-particle track.
  track.NewTrack(0, y0-0.00001, 0, 0, 0, -1, 0);
  // Retrieve the clusters along the track.
  for (const auto &cluster : track.GetClusters()) {
    // Loop over the electrons in the cluster.
    for (const auto &electron : cluster.electrons) {
      // Simulate an avalanche (until the set time window).
      aval.AvalancheElectron(electron.x, electron.y, electron.z, electron.t,
                             0., 0., 0., 0.);
      // Transfer the avalanche electrons to the grid.
      //avalgrid.AddElectrons(&aval);
    }
  }
  driftView.Plot();
  // Start grid based avalanche calculation starting from where the microsocpic
  // calculations stopped.
  //LOG("Switching to grid based method.");
  //avalgrid.StartGridAvalanche();
  // Stop timer.
  double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

  LOG("Script: " << "Electrons have drifted. It took " << duration
                 << "s to run.");

  if (plotSignal) {
    // Plot the induced current.
    signalView->PlotSignal(label);
    signalView->EnableLegend();
    cSignal->Update();
    gSystem->ProcessEvents();
    sensor.ExportSignal(label, "Signal");
    // Plot the induced charge.
    sensor.IntegrateSignal(label);
    chargeView->PlotSignal(label);
    cCharge->Update();
    gSystem->ProcessEvents();
    // Export induced current data as an csv file.
    sensor.ExportSignal(label, "Charge");
  }
  */
  /*
  LOG("Script: Total induced charge = " << sensor.GetTotalInducedCharge(label)
                                        << " [fC].");
  LOG("Field Area x = "<< -0.5 <<" to "<<0.5);
  LOG("Field Area y = "<< 0<<" to "<<dTotal);
  */
  // Export and plot Weighting Fields
  
/*
  if (plotField) {
    cField = new TCanvas("cField", "", 600, 600);
    fieldView = new ViewField(&sensor);
    
    fieldView->SetPlane(0,-y0,0,0,-0.001,0);
    fieldView->SetArea(-0.5,0,-0.5,0.5,dTotal,0.5);
    fieldView->SetVoltageRange(-160., 160.);
    
    fieldView->SetCanvas(cField);
    LOG("Calculate Strip Weighting Field - on the GRID");
    rpc->SetWeightingPotentialGrid(-0.5, 0.5, 1, 0, dTotal, 100, -0.5, 0.5, 1, label);
    LOG("Plot Strip Weighting Field - on the GRID");
    cField->SetLeftMargin(0.16);
    fieldView->PlotProfileWeightingField(label,0., 0., 0., 0., dTotal, 0.,"v",true);
    //fieldView->PlotWeightingField(label,"emag");
    //fieldView->PlotContour();  
  }
*/
  std::cin.get();
  LOG("End of Program");
  app.Run();
}
