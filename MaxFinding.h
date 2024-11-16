#pragma once

#include <algorithm>
#include <numeric>

#include <algorithms/algorithm.h>
#include <DD4hep/BitFieldCoder.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DDRec/Surface.h>
#include <DDRec/SurfaceManager.h>

#include <spdlog/spdlog.h>

// Event Model related classes
#include <edm4eic/CalorimeterHitCollection.h>
#include <edm4eic/ProtoClusterCollection.h>
#include <edm4hep/utils/vector_utils.h>

#include "algorithms/interfaces/WithPodConfig.h"
#include "MaxFindingConfig.h"

namespace eicrecon {

    using MaxFindingAlgorithm = algorithms::Algorithm<
        algorithms::Input<edm4eic::CalorimeterHitCollection>,
        algorithms::Output<edm4eic::CalorimeterHitCollection>
    >;

    class MaxFinding : public MaxFindingAlgorithm, 
			    public WithPodConfig<MaxFindingConfig> {
    public:
        MaxFinding(std::string_view name)
            : MaxFindingAlgorithm{name, {"inputHitCollection"}, {"outputMaximaCollection"},
                                   "Algorithm to find local energy maxima in calorimeter hits."} {}

    private:
        // Unitless inputs
	std::array<double, 2> localDistXY{0, 0};
	std::array<double, 2> layerDistXY{0, 0};
	double minHitEnergy{0};
	double peakEnergyThreshold{0};
	double neighborEnergyThreshold{0};
	int minNumofNeighbors{0};

    public:
        void init() {
        	if (m_cfg.localDistXY.size() != 2) {
            		error( "Expected 2 values (x_dist, y_dist) for localDistXY");
            		return;
        	}

		debug("Config - localDistXY: {}, {}", m_cfg.localDistXY[0], m_cfg.localDistXY[1]);
    		debug("Config - layerDistXY: {}, {}", m_cfg.layerDistXY[0], m_cfg.layerDistXY[1]);
    		debug("Config - minHitEnergy: {}", m_cfg.minHitEnergy);
    		debug("Config - peakEnergyThreshold: {}", m_cfg.peakEnergyThreshold);
    		debug("Config - neighborEnergyThreshold: {}", m_cfg.neighborEnergyThreshold);
    		debug("Config - minNumofNeighbors: {}", m_cfg.minNumofNeighbors);

        	// using juggler internal units (dd4hep:GeV, dd4hep::mm, dd4hep::ns, dd4hep::rad)
        	localDistXY[0] = m_cfg.localDistXY[0] / dd4hep::mm;
        	localDistXY[1] = m_cfg.localDistXY[1] / dd4hep::mm;
        	layerDistXY[0] = m_cfg.layerDistXY[0] / dd4hep::mm;
        	layerDistXY[1] = m_cfg.layerDistXY[1] / dd4hep::mm;
        	minHitEnergy = m_cfg.minHitEnergy / dd4hep::GeV;
		peakEnergyThreshold = m_cfg.peakEnergyThreshold / dd4hep::GeV;
        	neighborEnergyThreshold = m_cfg.neighborEnergyThreshold / dd4hep::GeV;
        	minNumofNeighbors = m_cfg.minNumofNeighbors;		
	
	}

	void process(const Input& input, const Output& output) const final {
    		const auto [hits] = input;  // Assuming 'hits' is a pointer to a CalorimeterHitCollection
    		auto& proto = std::get<0>(output);  // Access the CalorimeterHitCollection inside the tuple

    		// Iterate through all hits to find local maxima
    		for (const auto& hit : *hits) {
        		if(hit.getEnergy() <= minHitEnergy){continue;}
			if (is_local_maximum(hit, *hits)) {
				//std::cerr<<"!!!!!!!!!!!!!!!!Found peak hit!!!!!!!!!!!!"<<std::endl;
				auto newHit = proto->create();  // Create a new hit in the output collection
            			newHit.setEnergy(hit.getEnergy());   // Set energy from the original hit
            			newHit.setPosition(hit.getPosition()); // Set position from the original hit
            			newHit.setTime(hit.getTime());  // Set time from the original hit (if needed)
        		}
    		}

	}




    private:
	bool is_neighbor(const edm4eic::CalorimeterHit& h1, const edm4eic::CalorimeterHit& h2) const {
        	//if (h1.getEnergy() < minHitEnergy || h2.getEnergy() < minHitEnergy) return false;
		
		// layer check
        	int ldiff = std::abs(h1.getLayer() - h2.getLayer());

		// same layer, check local positions
        	if (ldiff == 0) {
			//std::cerr<<"Local distance in x: "<<std::to_string(std::abs(h1.getLocal().x - h2.getLocal().x)) << " local distance in y: " << std::to_string(std::abs(h1.getLocal().y - h2.getLocal().y)) <<std::endl;
			return (std::abs(h1.getLocal().x - h2.getLocal().x) <= localDistXY[0]) && (std::abs(h1.getLocal().y - h2.getLocal().y) <= localDistXY[1]);
        	} else if (ldiff <= m_cfg.neighborLayersRange) {
            		//std::cerr<<"global distance in x: "<<std::to_string(std::abs(h1.getPosition().x - h2.getPosition().x)) << " global distance in y: " << std::to_string(std::abs(h1.getPosition().y - h2.getPosition().y)) <<std::endl;
			return (std::abs(h1.getPosition().x - h2.getPosition().x) <= layerDistXY[0]) && (std::abs(h1.getPosition().y - h2.getPosition().y) <= layerDistXY[1]);
          	}
        	
        	// not in adjacent layers
        	return false;
    	}

	edm4eic::CalorimeterHitCollection getNeighbors(const edm4eic::CalorimeterHit& targetHit, const edm4eic::CalorimeterHitCollection& allHits) const {
    
    		edm4eic::CalorimeterHitCollection neighbors;
    
		for (const auto& hit : allHits) {
   			if (!(hit.getEnergy() == targetHit.getEnergy() && hit.getPosition() == targetHit.getPosition() && hit.getTime() == targetHit.getTime())) {
    				// Check if hit qualifies as a neighbor
    				if (hit.getEnergy() >= minHitEnergy && is_neighbor(targetHit, hit)) {
        				auto newHit = neighbors.create(); // Create a new hit in the new collection
        				newHit.setEnergy(hit.getEnergy());
        				newHit.setPosition(hit.getPosition());
        				newHit.setTime(hit.getTime());
    				}
			}
		}

    		return neighbors;
	}


	// Function to calculate the total energy of neighbors
	double calculateTotalNeighborEnergy(const edm4eic::CalorimeterHit& hit, const edm4eic::CalorimeterHitCollection& neighbors) const {
    		// Sum the energy of all neighboring hits
    		double totalNeighborEnergy = std::accumulate(neighbors.begin(), neighbors.end(), 0.0, [](double sum, const edm4eic::CalorimeterHit& neighbor) {return sum + neighbor.getEnergy();});

    		return totalNeighborEnergy;
	}


	// Function to pull total number of neighbors
	int totalNumberofNeighbors(const edm4eic::CalorimeterHit& hit, const edm4eic::CalorimeterHitCollection& allHits) const {
		auto neighbors = getNeighbors(hit, allHits);
		return neighbors.size();
	}

	bool is_local_maximum(const edm4eic::CalorimeterHit& hit, const edm4eic::CalorimeterHitCollection& hits) const {
		auto neighboringHits = getNeighbors(hit, hits);
		int totalNeighbors = neighboringHits.size();
    		double neighborEnergy = calculateTotalNeighborEnergy(hit, neighboringHits);
		
		//std::cerr<<"Total neighbors: " << std::to_string(totalNeighbors) <<" with energy: "<<std::to_string(neighborEnergy)<<std::endl;
		//std::cerr<<"Peak energy of: " << std::to_string(hit.getEnergy()) <<std::endl;

    		if (totalNeighbors < minNumofNeighbors || neighborEnergy < neighborEnergyThreshold) {
			//std::cerr<<"Failed neighbor num check"<<std::endl;
			return false;
    		}

    		if (hit.getEnergy() < peakEnergyThreshold) {
			//std::cerr<<"Failed peak energy threshold"<<std::endl;
        		return false;
    		}

    		for (const auto& other : neighboringHits) {
        		if (hit.getEnergy() <= other.getEnergy()) {
				//std::cerr<<"Is not local maximum" <<std::endl;
				//std::cerr<<"Num of neighbors: " << std::to_string(totalNeighbors) << std::endl;
				//std::cerr<<"Energy of this hit: " << std::to_string(hit.getEnergy()) << " energy of the other hit: " << std::to_string(other.getEnergy()) << std::endl;
				return false;
        		}
    		}
    		return true;
	}


    };

} // namespace eicrecon

