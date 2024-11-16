#pragma once

#include "algorithms/calorimetry/MaxFinding.h"
#include "extensions/jana/JOmniFactory.h"
#include "services/algorithms_init/AlgorithmsInit_service.h"

namespace eicrecon {

class PeakFinding_factory : public JOmniFactory<PeakFinding_factory, MaxFindingConfig> {

public:
    	using AlgoT = eicrecon::MaxFinding;
private:
	std::unique_ptr<AlgoT> m_algo;

    	PodioInput<edm4eic::CalorimeterHit> m_hits_input {this};
    	PodioOutput<edm4eic::CalorimeterHit> m_protos_output {this};

    	ParameterRef<std::vector<double>> m_ldxy {this, "localDistXY", config().localDistXY};
    	ParameterRef<std::vector<double>> m_ldxy_adjacent {this, "layerDistXY", config().layerDistXY};
   	ParameterRef<double> m_mhe {this, "minHitEnergy", config().minHitEnergy};
	ParameterRef<double> m_pet {this, "peakEnergyThreshold", config().peakEnergyThreshold};
	ParameterRef<double> m_net {this, "neighborEnergyThreshold", config().neighborEnergyThreshold};
	ParameterRef<int> m_mnn {this, "minNumofNeighbors", config().minNumofNeighbors};

    	Service<AlgorithmsInit_service> m_algorithmsInit {this};

public:
    	void Configure() {
		//std::cerr<< "Inside config"<<std::endl;
        	spdlog::set_level(spdlog::level::debug);
		
		if (logger() == nullptr) {
        		std::cerr << "Error: logger() returned nullptr." << std::endl;
        		return;  // Exit Configure to avoid using a null logger
    		}
    
    		//spdlog::debug("Logger is not nullptr; proceeding with configuration.");
		
		m_algo = std::make_unique<AlgoT>(GetPrefix());
        	m_algo->level(static_cast<algorithms::LogLevel>(logger()->level()));
        	m_algo->applyConfig(config());
        	m_algo->init();
    	}
		
    	void ChangeRun(int64_t run_number) {
    	}

    	void Process(int64_t run_number, uint64_t event_number) {
		 //m_algo->process({m_hits_input()}, {m_protos_output().get()});
		//std::cerr << "Processing event " << event_number << "..." << std::endl;

    		// Before processing, check if the input hits collection is empty or not
    		auto hits = m_hits_input();
    		//std::cerr << "Number of hits: " << hits->size() << std::endl;

    		// Process the hits and fill the protos output
    		m_algo->process({hits}, {m_protos_output().get()});

    		// After processing, check the output collection
    		auto& protos = *m_protos_output();
    		//std::cerr << "Number of protos output: " << protos.size() << std::endl;
	}
};

} // namespace eicrecon
