#pragma once

#include <string>
#include <iostream>
#include <vector>

namespace eicrecon {

  struct MaxFindingConfig {

    // maximum distance of local (x, y) to be considered as neighbors at the same layer
    std::vector<double> localDistXY = {1.0 * dd4hep::mm, 1.0 * dd4hep::mm};
    // maximum distance of global (x, y) to be considered as neighbors at different layers
    std::vector<double> layerDistXY = {1.0 * dd4hep::mm, 1.0 * dd4hep::mm};

    // Number of neighboring layers to consider nieghboring for a hit
    int neighborLayersRange = 1;

    // minimum hit energy to participate clustering
    double minHitEnergy = .05 * dd4hep::MeV;
    //minimum energy for maximum
    double peakEnergyThreshold = 5 * dd4hep::MeV;
    // minimum cluster energy
    double neighborEnergyThreshold = 10 * dd4hep::MeV;
    // minimum number of neighboring hits
    int minNumofNeighbors = 4;

  };  
} // namespace eicrecon
