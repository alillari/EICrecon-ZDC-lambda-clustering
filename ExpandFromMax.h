#pragma once

#include <algorithm>
#include <numeric>

#include <optional>
#include <unordered_map>

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

    using ExpandFromMax = algorithms::Algorithm<
        algorithms::Input<edm4eic::CalorimeterHitCollection, edm4eic::CalorimeterHitCollection>, // Two inputs
        algorithms::Output<edm4eic::ProtoClusterCollection>
    >;

    class ExpandFromMax : public ExpandFromMaxAlgorithm, 
			    public WithPodConfig<ExpandFromMaxConfig> {
    public:
        ExpandFromMax(std::string_view name)
            : ExpandFromMaxAlgorithm{name, {"inputHitCollection", "inputMaxHitCollection"}, {"outputProtoClusters"},
                                   "Algorithm to create local clusters around local energy maximums."} {}

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

	}




    private:
	bool is_neighbor(const edm4eic::CalorimeterHit& h1, const edm4eic::CalorimeterHit& h2) const {
		// layer check
        	int ldiff = std::abs(h1.getLayer() - h2.getLayer());

		// same layer, check local positions
        	if (ldiff == 0) {
			return (std::abs(h1.getLocal().x - h2.getLocal().x) <= localDistXY[0]) && (std::abs(h1.getLocal().y - h2.getLocal().y) <= localDistXY[1]);
        	} else if (ldiff <= m_cfg.neighborLayersRange) {
			return (std::abs(h1.getPosition().x - h2.getPosition().x) <= layerDistXY[0]) && (std::abs(h1.getPosition().y - h2.getPosition().y) <= layerDistXY[1]);
          	}
        	
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

	double EuclideanDistance(const edm4eic::CalorimeterHit& hit1, const edm4eic::CalorimeterHit& hit2) {
    		return std::sqrt(std::pow(hit1.getPosition().x - hit2.getPosition().x, 2) + std::pow(hit1.getPosition().y - hit2.getPosition().y, 2) + std::pow(hit1.getPosition().z - hit2.getPosition().z, 2));
	}	

	std::vector<size_t> getEucldNeighbors(const edm4eic::CalorimeterHit& targetHit, const edm4eic::CalorimeterHitCollection& allHits, double minDist) {
		std::vector<size_t> neighbors;

    		for (size_t int i = 0; i < allHits.size(); i++) {
        		if (!(allHits[i].getEnergy() == targetHit.getEnergy() && allHits[i].getPosition() == targetHit.getPosition() && allHits[i].getTime() == targetHit.getTime())){
				if (EuclideanDistance(targetHit, allHits[i]) <= minDist) {
        				neighbors.push_back(i);
				}
			}
    		}

    		return neighbors;
	}

	std::unordered_map<int, edm4eic::CalorimeterHitCollection> DBSCAN(edm4eic::CalorimeterHitCollection& hits, double minDist, int minPts) {
    		int cluster_id = 1; // Start cluster IDs from 1 (0 is reserved for noise)
    		std::vector<int> cluster_ids(hits.size(), -1); // -1: unclassified, 0: noise
    		std::unordered_map<int, edm4eic::CalorimeterHitCollection> clusters;

    		for (size_t i = 0; i < hits.size(); ++i) {
        		if (cluster_ids[i] != -1) continue; // Skip already classified hits

        		// Find neighbors of the current hit
        		auto neighbors = getEucldNeighbors(hits, i, epsilon);

        		if (neighbors.size() < minPts) {
            			cluster_ids[i] = 0; // Mark as noise
        		} else {
            			// Create a new cluster
            			edm4eic::CalorimeterHitCollection cluster_hits;
            			ExpandCluster(hits, i, neighbors, cluster_id, epsilon, minPts, cluster_ids, cluster_hits);
            			clusters[cluster_id] = cluster_hits; // Add cluster to map
            			++cluster_id;
        		}
    		}

    		return clusters;	
	}
	
	void ExpandCluster(const edm4eic::CalorimeterHitCollection& hits, size_t index, std::vector<size_t> neighbors, int cluster_id, double epsilon, int minPts, std::vector<int>& cluster_ids) {
		cluster_ids[index] = cluster_id;
    		cluster_hits.push_back(hits[index]);

    		for (size_t i = 0; i < neighbors.size(); ++i) {
        		size_t neighbor_index = neighbors[i];

        		if (cluster_ids[neighbor_index] == 0) {
            			cluster_ids[neighbor_index] = cluster_id; // Convert noise to cluster member
        		}

        		if (cluster_ids[neighbor_index] != -1){continue;} // Skip already classified hits

        		cluster_ids[neighbor_index] = cluster_id;
        		cluster_hits.push_back(hits[neighbor_index]);

        		// Recursively find and add neighbors
        		auto new_neighbors = GetNeighbors(hits, neighbor_index, epsilon);
        		if (new_neighbors.size() >= minPts) {
            			neighbors.insert(neighbors.end(), new_neighbors.begin(), new_neighbors.end());
        		}
    		}
	}

	double getClusterEnergy(const edm4eic::CalorimeterHitCollection& cluster){
		double cluster_energy = 0;

		for(const auto& hit: cluster){
			cluster_energy = cluster_energy + hit.getEnergy();
		}

		return cluster_energy;
	}

	std::optional<std::unordered_map<int, edm4eic::CalorimeterHitCollection>> sortMaxClusters(std::unordered_map<int, edm4eic::CalorimeterHitCollection> max_clusters){
		std::unordered_map<int, edm4eic::CalorimeterHitCollection> copied_clusters(max_clusters)
		
		if(max_clusters.size() != 3){
			return std::nullptr;
		}else{
			double cluster1_energy = getClusterEnergy(max_clusters[1]);
			double cluster2_energy = getClusterEnergy(max_clusters[2]);
			double cluster3_energy = getClusterEnergy(max_clusters[3]);

			if(cluster1_energy > cluster2_energy && cluster1_energy > cluster3_energy){
				copied_clusters.erase(1);
			} else if(cluster2_energy > cluster1_energy && cluster2_energy > cluster3_energy){
				copied_clusters.erase(2);
			} else if(cluster3_energy > cluster1_energy && cluster3_energy > cluster2_energy){
				copied_clusters.erase(3);
			}
			return copied_clusters;
		}
	}


	

    };

} // namespace eicrecon

