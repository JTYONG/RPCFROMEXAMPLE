#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <ctime>
#include <iostream>
#include <numeric>
#include <vector>

#include "Garfield/AvalancheGrid.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentParallelPlate.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"

#define LOG(x) std::cout << x << std::endl

using namespace Garfield;

int main(int argc, char *argv[]) {
  // Creates ROOT Application Environtment
  TApplication app("app", &argc, argv);

  const bool debug = true;
  constexpr bool plotSignal = true;
  constexpr bool plotField = true;
  
  // Describe the geometry of the RPC.
  // Relative permitivity of the layers
  const double epMylar = 3.1;     // [1]
  const double epGlaverbel = 8.;  // [1]
  const double epWindow = 6.;     // [1]
  const double epGas = 1.;        // [1]
  std::vector<double> eps = {epMylar, epWindow,    epGas,  epGlaverbel,
                             epGas,   epGlaverbel, epGas,  epGlaverbel,
                             epGas,   epGlaverbel, epGas,  epGlaverbel,
                             epGas,   epWindow,    epMylar};

  // Thickness of the layers
  const double dMylar = 0.035;     // [cm]
  const double dGlaverbel = 0.07;  // [cm]
  const double dWindow = 0.12;     // [cm]
  const double dGas = 0.025;       // [cm]
  std::vector<double> thickness = {dMylar, dWindow,    dGas,  dGlaverbel,
                                   dGas,   dGlaverbel, dGas,  dGlaverbel,
                                   dGas,   dGlaverbel, dGas,  dGlaverbel,
                                   dGas,   dWindow,    dMylar};
  const double dTotal = std::accumulate(thickness.begin(), thickness.end(), 0.);

  std::cout<<"Total Thickness: "<<2*dMylar + 2*dWindow + 6*dGas +5*dGlaverbel<<std::endl;

  // Applied potential
  const double voltage = 15e3;  // [V] = 15 kV

  ComponentParallelPlate *rpc = new ComponentParallelPlate();
  rpc->Setup(eps.size(), eps, thickness, voltage);

  // Adding a readout structure.
  const std::string label = "ReadoutPlane";
  rpc->AddPlane(label);

  // Set up the gas (C2H2F4/iC4H10/SF6 90/5/5).
  MediumMagboltz gas;
  gas.LoadGasFile("c2h2f4_ic4h10_sf6.gas");
  gas.Initialise(true);

  // Set the drift medium.
  rpc->SetMedium(&gas);

  // Create the sensor.
  Sensor sensor(rpc);
  sensor.AddElectrode(rpc, label);

  // Set the time bins.
  const std::size_t nTimeBins = 200;
  const double tmin = 0.;
  const double tmax = 4;
  const double tstep = (tmax - tmin) / nTimeBins;
  sensor.SetTimeWindow(tmin, tstep, nTimeBins);

  // Use microscopic tracking for initial stage of the avalanche.
  AvalancheMicroscopic aval(&sensor);
  // Set time until which the calculations will be done microscopically.
  const double tMaxWindow = 0.1;
  aval.SetTimeWindow(0., tMaxWindow);

  // Use a grid-based method for simulating the avalanche growth 
  // after the initial stage.
  AvalancheGrid avalgrid(&sensor);

  const int nX = 5;
  const int nZ = 5;
  const double dY = 1.e-4;
  const int nY = int(dTotal / dY);
  avalgrid.SetGrid(-0.05, 0.05, nX, 0., dTotal, nY, -0.05, 0.05, nZ);

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
  track.SetMomentum(1.e11);
  track.CrossInactiveMedia(true);

  std::clock_t start = std::clock();

  auto cD = new TCanvas("cD","Drift",600,600);
  ViewDrift driftView;
  driftView.SetCanvas(cD);
  aval.EnablePlotting(&driftView, 100);
  track.EnablePlotting(&driftView);

  // Start the track in the first gas layer.
  const double y0 = dTotal - dMylar - dWindow;
  // Simulate a charged-particle track.
  track.NewTrack(0, y0, 0, 0, 0, -1, 0);
  // Retrieve the clusters along the track.
  for (const auto &cluster : track.GetClusters()) {
    // Loop over the electrons in the cluster.
    for (const auto &electron : cluster.electrons) {
      // Simulate an avalanche (until the set time window).
      aval.AvalancheElectron(electron.x, electron.y, electron.z, electron.t,
                             0.1, 0., 0., 0.);
      // Transfer the avalanche electrons to the grid.
      avalgrid.AddElectrons(&aval);
    }
  }
  driftView.Plot();
  // Start grid based avalanche calculation starting from where the microsocpic
  // calculations stopped.
  LOG("Switching to grid based method.");
  avalgrid.StartGridAvalanche();
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
  /*
  LOG("Script: Total induced charge = " << sensor.GetTotalInducedCharge(label)
                                        << " [fC].");
  LOG("Field Area x = "<< -0.5 <<" to "<<0.5);
  LOG("Field Area y = "<< 0<<" to "<<dTotal);
  */
  // Export and plot Weighting Fields
  ViewField *fieldView = nullptr;
  TCanvas *cField = nullptr;
  
  if (plotField) {
    cField = new TCanvas("cField", "", 600, 600);
    fieldView = new ViewField(&sensor);
    /*
    fieldView->SetPlane(0,-y0,0,0,-0.001,0);
    fieldView->SetArea(-0.5,0,-0.5,0.5,dTotal,0.5);
    fieldView->SetVoltageRange(-160., 160.);
    */
    fieldView->SetCanvas(cField);
    LOG("Calculate Strip Weighting Field - on the GRID");
    rpc->SetWeightingPotentialGrid(-0.5, 0.5, 1, 0, dTotal, 100, -0.5, 0.5, 1, label);
    LOG("Plot Strip Weighting Field - on the GRID");
    cField->SetLeftMargin(0.16);
    fieldView->PlotProfileWeightingField(label,0., 0., 0., 0., dTotal, 0.,"v",true);
    //fieldView->PlotWeightingField(label,"emag");
    //fieldView->PlotContour();  
  }

  LOG("End of Program");
  app.Run(true);
}
