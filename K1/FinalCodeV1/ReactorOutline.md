# Bioreactor Writeup Outline
* Brief introduction
m* CHO metabolic model
	* Reaction set
	* Kinetic rates
	* MFB to get all rates
	 	* Quadratic programming in MATLAB, tolerance etc
	* RK4 to model changing concentrations over time
	* Explanation of R, redox parameter
* Issues with model
	* Poor numerical method
	 	* R not matching with lower step size
ed	* Refitting the parameters
ed	 	* Our multiple objective functions
ed	 	* Our multiple fitting methods (particle swarm, simulated annealing, and pattern search)
ed	 	* Final methods: minimize SSR with particle swarm to get approximate global min, then pattern search
ed		* Comment on agreement & error
ed* New parameters, batch output
ed	* Overall reaction, production, graphs, etc
* Generalizations of the model
m	* Adaptation for perfusion
		* Overall reaction, production, graphs, etc
ed	* Adaptation for the seed train
ed		* Overall reaction, production, graphs, etc 
m* O2/CO2 systems
	* Reactor specifications
	 	* Sparger type, diameters, heights, impellers, etc
	* Model for bubble properties
	 	* KLa, eps, P/V etc
	* Control system
	 	* Return with small time step
	 	* Describe PI controllers, multiple gas feeds
			* valve specifications
	 	* Outputs, gas flow rates required
		* Impeller control system, eddy size
ed* Heat transfer
ed	* Heat transfer model
ed		* bubble transfer, jacker transfer, cell heat generation
ed	* Model parameters for each vessel
ed	* outputs (utility requirements, min/max temperatures in each vessel)
* Summary and conclusions (brief)
* Summary of Assumptions
