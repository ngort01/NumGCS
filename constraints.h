#pragma once
class constraints
{
public:
	constraints();
	~constraints();
	typedef double(*f)(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
		std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices);

	static double ec(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double dr(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double dc(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double tangent_internal(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double non_tangent_internal(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double collinear(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double perpendicular(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double parallel(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double coincident(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double angle(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double diff(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double angle90cc(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double angle_eq(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double leq(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double lt(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double midpoint(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	static double arc_length(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	//static double distance_pl(const std::vector<double>&, std::vector<double>&, std::vector<int>&, std::vector<bool>&, std::vector<double>&, std::map<int, int>&, std::map<int, int>&);
	
	f get_obj_func(const std::string& name) {
		if (funcMap[name]) {
			return funcMap[name];
		}
		else {
			std::cerr << "Constraint not found!" << std::endl;
		}
	}

private:
	std::map<std::string, f> funcMap;
	
};


