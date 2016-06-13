#include "stdafx.h"
#include "lib.h"
#include"constraints.h"

static constraints constraint_ref;
static lib* _lib;

lib::lib()
{
	maxEval = 50000;
	evalCount = 0;
	threshold = 0.01;
	main_alg = "G_MLSL"; // main algorithm, global or local
	//main_alg = "LD_LBFGS"; // main algorithm, global or local
	local_alg = "LN_BOBYQA"; // local algorithm, if main is global
	//local_alg = "LN_COBYLA"; // local algorithm, if main is global
	//local_alg = "LD_LBFGS"; // local algorithm, if main is global

	// local gradient based
	algMap["LD_LBFGS"] = nlopt::LD_LBFGS;
	algMap["LD_MMA"] = nlopt::LD_MMA;
	algMap["LD_CCSAQ"] = nlopt::LD_CCSAQ;
	algMap["LD_SLSQP"] = nlopt::LD_SLSQP;
	algMap["LD_TNEWTON"] = nlopt::LD_TNEWTON;
	algMap["LD_VAR1"] = nlopt::LD_VAR1;
	algMap["LD_VAR2"] = nlopt::LD_VAR2;

	// local derivative-free
	algMap["LN_BOBYQA"] = nlopt::LN_BOBYQA;
	algMap["LN_COBYLA"] = nlopt::LN_COBYLA;
	algMap["LN_SBPLX"] = nlopt::LN_SBPLX;

	// global
	algMap["G_MLSL"] = nlopt::G_MLSL;
	algMap["G_MLSL_LDS"] = nlopt::G_MLSL_LDS;
	algMap["GN_DIRECT_L"] = nlopt::GN_DIRECT_L;
	algMap["GN_DIRECT_L_RAND"] = nlopt::GN_DIRECT_L_RAND;
	algMap["GN_CRS2_LM"] = nlopt::GN_CRS2_LM;
	algMap["GD_STOGO"] = nlopt::GD_STOGO;
	algMap["GD_STOGO_RAND"] = nlopt::GD_STOGO_RAND;
	algMap["GN_ISRES"] = nlopt::GN_ISRES;
	algMap["GN_ESCH"] = nlopt::GN_ESCH;

	algMap["AUGLAG"] = nlopt::AUGLAG;

	_lib = this;
}


lib::~lib()
{
}

////////////////////////////////////////////////////////////////
//////////////////// add parameters ////////////////////////////
////////////////////////////////////////////////////////////////
int lib::addParameter(double v, double lb, double ub, bool f)
{
	all_params.push_back(v);

	lBounds.push_back(lb);
	uBounds.push_back(ub);
	fixed.push_back(f);
	return all_params.size() - 1;
}

std::vector<int> lib::addParameters(std::vector<double>& vals, std::vector<double>& lb, std::vector<double>& ub, std::vector<bool>& f)
{
	std::vector<int> indices;
	for (std::vector<double>::size_type i = 0; i != vals.size(); i++) {
		indices.push_back(addParameter(vals[i], lb[i], ub[i], f[i]));
	}
	return indices;
}


////////////////////////////////////////////////////////////////
//////////////////// get parameters ////////////////////////////
////////////////////////////////////////////////////////////////
double lib::getParameter(int index)
{
	return all_params[index];
}

std::vector<double> lib::getParameters(std::vector<int> indices)
{
	int size = indices.size();
	std::vector<double> p(size);
	for (int i = 0; i != size; i++) {
		p[i] = all_params[i];
	}
	return p;
}

std::vector<double> lib::getAllParameters()
{
	return all_params;
}


////////////////////////////////////////////////////////////////
//////////////////// add constraints ///////////////////////////
////////////////////////////////////////////////////////////////
void lib::addConstraint(std::string id, std::vector<int>& indices) {
	cons[id].push_back(indices);
}



////////////////////////////////////////////////////////////////
//////////////////// set algorithms ////////////////////////////
////////////////////////////////////////////////////////////////
void lib::setAlgorithms(std::string main)
{
	main_alg = main;
}
void lib::setAlgorithms(std::string main, std::string local)
{
	main_alg = main;
	local_alg = local;
}

////////////////////////////////////////////////////////////////
////////////////////// set max eval ////////////////////////////
////////////////////////////////////////////////////////////////
void lib::setMaxEval(int n) {
	maxEval = n;
}

double lib::getMinVal()
{
	return minf;
}

void lib::setThreshold(double t)
{
	threshold = t;
}

int lib::getEvalCount()
{
	return evalCount;
}


////////////////////////////////////////////////////////////////
////////////////////////// reset ///////////////////////////////
////////////////////////////////////////////////////////////////
void lib::clearSolver()
{

	//std::cout << "Clearing solver ..." << std::endl;
	all_params.clear();

	//var_params.clear();
	var_indices.clear();

	fixed_params.clear();
	fixed_indices.clear();

	fixed.clear();

	lBounds.clear();
	uBounds.clear();
	cons.clear();

	//std::cout << "Clearing complete."<< std::endl;

}


////////////////////////////////////////////////////////////////
//////////////////// objective function ////////////////////////
////////////////////////////////////////////////////////////////


double obj_func(const std::vector<double> &x, std::vector<double> &grad, void* data) {
	++_lib->evalCount;
	double f = 0.0;

	for (std::map<std::string, std::vector<std::vector<int> > >::iterator obj = _lib->cons.begin();
	obj != _lib->cons.end();
		obj++) {

		std::string key = (*obj).first;

		//for (auto& indices : obj.second) {
		for (std::vector<std::vector<int> >::iterator indices = (*obj).second.begin();
		indices != (*obj).second.end();
			indices++) {
			f += constraint_ref.get_obj_func(key)(x, grad, *indices, _lib->fixed, _lib->fixed_params, _lib->fixed_indices, _lib->var_indices);
		}
	}

	return f;
}


////////////////////////////////////////////////////////////////
/////////////////////////// solve //////////////////////////////
////////////////////////////////////////////////////////////////
bool lib::solve()
{

	//- build vector of variable parameters
	std::vector<double> vars;
	std::vector<double> lb;
	std::vector<double> ub;
	for (int i = 0; i < fixed.size(); i++)
	{
		if (!fixed[i])
		{
			vars.push_back(all_params[i]);
			lb.push_back(lBounds[i]);
			ub.push_back(uBounds[i]);
			var_indices[i] = vars.size() - 1;
		}
		else
		{
			fixed_params.push_back(all_params[i]);
			fixed_indices[i] = fixed_params.size() - 1;
		}
	}


	int dim = vars.size();

	//std::cout << "DIM = " << dim << std::endl;

	if (dim == 0)
	{   //- no variable parameters, nothing to solve
		return true;
	}


	if (dim != lb.size() || dim != ub.size())
	{
		std::cerr << "NLOPT ERROR: number of variable parameters (" << dim << ") does not equal number of bounded variable parameters (lower bound count = " << lb.size() << ", and upper bound count = " << ub.size() << "). Aborting solve." << std::endl;
		return false;

	}

	/* algorithm and dimensionality */
	//nlopt::opt* opt = new nlopt::opt(algMap[main_alg], dim);
	nlopt::opt opt(algMap[main_alg], dim);
	nlopt::opt local_opt(algMap[local_alg], dim);


	std::cout << "Algorithms: " << main_alg << " - " << local_alg << std::endl;

	/*bounds*/
	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);

	/* objective function */
	opt.set_min_objective(obj_func, NULL);

	opt.set_xtol_rel(1e-8);
	opt.set_stopval(threshold);
	opt.set_maxeval(maxEval);

	opt.set_local_optimizer(local_opt);


	double minf = HUGE_VAL;
	nlopt::result result;

	//    NLOPT_SUCCESS = 1, /* generic success code */
	//    NLOPT_STOPVAL_REACHED = 2,
	//    NLOPT_FTOL_REACHED = 3,
	//    NLOPT_XTOL_REACHED = 4,
	//    NLOPT_MAXEVAL_REACHED = 5,
	//    NLOPT_MAXTIME_REACHED = 6

	try
	{
		result = opt.optimize(vars, minf);
	}
	catch (nlopt::roundoff_limited& e)
	{
		//- no more optimisation progress due to round off errors
		std::cout << "NLOPT FAIL (due to round off limits): " << e.what() << std::endl;
		return false;
	}
	catch (std::runtime_error& e)
	{
		//- failure
		std::cout << "NLOPT FAIL: " << e.what() << ", Min: " << opt.last_optimum_value() << std::endl;
		return false;
	}
	catch (std::exception& e)
	{
		std::cerr << "An NLOPT exception occurred: " << e.what() << std::endl;
		clearSolver();

		return false;
	}

	// code 5 returns true but the result may be wrong
	if (minf > threshold) {
		std::cout << "Min value too high! min: " << minf << " Threshold: " << threshold << std::endl;
		return false;
	}


	// on succes add optimized parameters to the all_params vector
	for (std::map<int, int>::iterator obj = var_indices.begin(); obj != var_indices.end(); obj++) {
		int key = (*obj).first;
		int var_idx = (*obj).second;
		all_params[key] = vars[var_idx];
	}

	std::cout << "Max eval: " << opt.get_maxeval() << std::endl;

	std::cout << "Min: " << minf << ", Result: " << result << std::endl;

	return true;

}
