#pragma once

class lib
{
public:
	lib();
	~lib();
	/*
	Adds a parameter(s) to the parameter vector.
		Input:
			double: parameter(s) 
			double: lower bound(s) 
			double: upper bound(s) 
			bool: fixed

		Returns:
			index/indices of the added parameter
	*/
	int addParameter(double, double, double, bool);
	std::vector<int> addParameters(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<bool>&);

	/*
	Gets a Parameter by index.
	*/
	double getParameter(int);
	std::vector<double> getParameters(std::vector<int>);
	std::vector<double> getAllParameters();

	/*
	Adds a constraint.
		Input:
			string: constraint id (e.g. "ec")
			vector<int>: indices of the parameters
	*/
	void addConstraint(std::string, std::vector<int>&);

	/*
	Set algorithms for solving
		Input:
			main algorithm: can be global or local. default is "G_MLSL"
			local algorithm: if the main algorithm is glaobal. default is "LN_BOBYQA"
	*/
	void setAlgorithms(std::string);
	void setAlgorithms(std::string, std::string);

	void setMaxEval(int); // set max iterations
	double getMinVal(); // get the min value
	void setThreshold(double);
	int getEvalCount();

	/*
	Clears params and so on
	*/
	void clearSolver();

	/*
	Solves the contraints.
		Returns:
			bool: true -> success 
	*/

	bool solve();


	std::map<std::string, nlopt::algorithm> algMap; // map of algorithms
	std::string main_alg; // main algorithm, global or local
	std::string local_alg; // local algorithm, if main is global

	std::vector<double> all_params; // all parameters
	std::vector<bool> fixed; // flag showing if param is fixed

	//std::vector<double> var_params; // variable parameters
	std::map<int, int > var_indices; // mapping index in all_params to index in var_params

	std::vector<double> fixed_params; // fixed parameters
	std::map<int, int> fixed_indices; 

	std::vector<double> lBounds; // lower bounds
	std::vector<double> uBounds; // upper bounds

	std::map<std::string, std::vector<std::vector<int> > > cons; // map of constraints
	int evalCount;

private:
	int maxEval; // max number of evaluations
	double threshold;
	double minf;
};

