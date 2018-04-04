# Bioreactor Writeup Outline
* Brief introduction
* CHO metabolic model
	* Reaction set
	* Kinetic rates
	* MFB to get all rates
	 	* Quadratic programming in MATLAB, tolerance etc
	* RK4 to model changing concentrations over time
	* Explanation of R, redox parameter
* Issues with model
	* Poor numerical method
	 	* R not matching with lower step size
	* Refitting the parameters
	 	* Our multiple objective functions
	 	* Our multiple fitting methods (particle swarm, simulated annealing, and pattern search)
	 	* Final methods: minimize SSR with particle swarm to get approximate global min, then pattern search
		* Comment on agreement & error
* New parameters, batch output
	* Overall reaction, production, graphs, etc
* Generalizations of the model
	* Adaptation for perfusion
		* Overall reaction, production, graphs, etc
	* Adaptation for the seed train
		* Overall reaction, production, graphs, etc 
* O2/CO2 systems
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
* Heat transfer
	* Heat transfer model
		* bubble transfer, jacker transfer, cell heat generation
	* Model parameters for each vessel
	* outputs (utility requirements, min/max temperatures in each vessel)
* Summary and conclusions (brief)
* Summary of Assumptions