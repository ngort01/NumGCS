#include "stdafx.h"
#include "constraints.h"


constraints::constraints()
{
	funcMap["ec"] = constraints::ec;
	funcMap["dc"] = constraints::dc;
	funcMap["dr"] = constraints::dr;
	funcMap["tangent_internal"] = constraints::tangent_internal;
	funcMap["non_tangent_internal"] = constraints::non_tangent_internal;
	funcMap["collinear"] = constraints::collinear;
	funcMap["perpendicular"] = constraints::perpendicular;
	funcMap["parallel"] = constraints::parallel;
	funcMap["angle"] = constraints::angle;
	funcMap["diff"] = constraints::diff;
	funcMap["coincident"] = constraints::coincident; 
	funcMap["angle90cc"] = constraints::angle90cc;
	funcMap["angle_eq"] = constraints::angle_eq;
    funcMap["leq"] = constraints::leq;
    funcMap["lt"] = constraints::lt;
    funcMap["midpoint"] = constraints::midpoint;
    funcMap["arc_length"] = constraints::arc_length;

}


constraints::~constraints()
{
}

/*
 distance_pl: distance from point p to line (p1,p2)
  ...work in progress
 */

//double constraints::distance_pl(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
//                        std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
//{
//    //- derivative (d1*d1 * ((x2-x1)^2 + (y2-y1)^2) - (-x0*(y2-y1)+y0*(x2-x1)+x1*y2-x2*y1)^2)^2
//    
//	int ixp = indices[0], iyp = indices[1];
//    int ix1 = indices[2], iy1 = indices[3];
//    int ix2 = indices[4], iy2 = indices[5];
//    int idist = indices[6];
//    
//	double xp = fixed[ixp] ? fixed_params[fixed_indices[ixp]] : x[var_indices[ixp]];
//	double yp = fixed[iyp] ? fixed_params[fixed_indices[iyp]] : x[var_indices[iyp]];
//	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];
//	double y1 = fixed[iy1] ? fixed_params[fixed_indices[iy1]] : x[var_indices[iy1]];
//	double x2 = fixed[ix2] ? fixed_params[fixed_indices[ix2]] : x[var_indices[ix2]];
//	double y2 = fixed[iy2] ? fixed_params[fixed_indices[iy2]] : x[var_indices[iy2]];
//	double dist = fixed[idist] ? fixed_params[fixed_indices[idist]] : x[var_indices[idist]];
//    
//    
//    double dx = x2 - x1;
//    double dy = y2 - y1;
//    double len2 = dx*dx + dy*dy;
//    double rawarea = (-xp*dy+yp*dx+x1*y2-x2*y1);
//    double area2 = rawarea * rawarea;
//    
//    
//	double fx = dist*dist * len2 - area2;
//    double fx2 = fx*fx;
//    
//	if (!grad.empty()) {
//		if (!fixed[ixp]) grad[var_indices[ixp]] = -4*(y1 - y2)*(x2*(yp - y1) + xp * (y1 - y2) + x1 * (-yp + y2)) * (-(x2 * (-yp + y1) + x1 * (yp - y2) + xp * (-y1 + y2))^2 + dist^2 ((-x1 + x2)^2 + (-y1 + y2)^2));
//		if (!fixed[iyp]) grad[var_indices[iyp]] = -4*(-x1 + x2)*(x2*(yp - y1) + xp * (y1 - y2) + x1 * (-yp + y2)) * (-(x2 * (-yp + y1) + x1 * (yp - y2) + xp (-y1 + y2))^2 + dist^2 ((-x1 + x2)^2 + (-y1 + y2)^2));
//	}
//    
//	return  fx2;
//}




/*
 arc_length: length of arc of circle with radius r spanning theta radians
 */

double constraints::arc_length(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
                             std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
    //- derivative (l - r * t)^2
    
	int ilen   = indices[0];
	int ir     = indices[1];
	int itheta = indices[2];
    
	double len   = fixed[ilen] ? fixed_params[fixed_indices[ilen]] : x[var_indices[ilen]];
	double r     = fixed[ir] ? fixed_params[fixed_indices[ir]] : x[var_indices[ir]];
	double theta = fixed[itheta] ? fixed_params[fixed_indices[itheta]] : x[var_indices[itheta]];
    
    double fx = len - r * theta;
    
	double fx2 = fx*fx;
    
	if (!grad.empty()) {
		if (!fixed[ir])     grad[var_indices[ir]]     = 2 * theta * (r * theta - len);
		if (!fixed[itheta]) grad[var_indices[itheta]] = 2 * r * (r * theta - len);
		if (!fixed[ilen])   grad[var_indices[ilen]]   = 2 * (len - r * theta);
	}
    
	return  fx2;
}




/*
 midpoint: p0 is midpoint of line (p1, p2)
 */

double constraints::midpoint(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
                        std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
    //- derivative (x0 - (x1 + x2) / 2)^2
    
	int ix0 = indices[0], iy0 = indices[1];
	int ix1 = indices[2], iy1 = indices[3];
	int ix2 = indices[4], iy2 = indices[5];
    
	double x0 = fixed[ix0] ? fixed_params[fixed_indices[ix0]] : x[var_indices[ix0]];
	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];
	double x2 = fixed[ix2] ? fixed_params[fixed_indices[ix2]] : x[var_indices[ix2]];
	double y0 = fixed[iy0] ? fixed_params[fixed_indices[iy0]] : x[var_indices[iy0]];
	double y1 = fixed[iy1] ? fixed_params[fixed_indices[iy1]] : x[var_indices[iy1]];
	double y2 = fixed[iy2] ? fixed_params[fixed_indices[iy2]] : x[var_indices[iy2]];
    
    double fx = x0 - (x1 + x2) / 2;
    double fy = y0 - (y1 + y2) / 2;
    
	double fx2 = fx*fx + fy*fy;
    
	if (!grad.empty()) {
		if (!fixed[ix0]) grad[var_indices[ix0]] = 2*x0 - x1 - x2;
		if (!fixed[ix1]) grad[var_indices[ix1]] = 0.5 * (-2 * x0 + x1 + x2);
		if (!fixed[ix2]) grad[var_indices[ix2]] = 0.5 * (-2 * x0 + x1 + x2);
		if (!fixed[iy0]) grad[var_indices[iy0]] = 2*y0 - y1 - y2;
		if (!fixed[iy1]) grad[var_indices[iy1]] = 0.5 * (-2 * y0 + y1 + y2);
		if (!fixed[iy2]) grad[var_indices[iy2]] = 0.5 * (-2 * y0 + y1 + y2);
	}
    
	return  fx2;
}


/*
 LEQ: x0 =< x1
 */

double constraints::leq(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
                        std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
    //- derivative (x1 - x0)^2
    
	int ix0 = indices[0];
	int ix1 = indices[1];
    
	double x0 = fixed[ix0] ? fixed_params[fixed_indices[ix0]] : x[var_indices[ix0]];
	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];
    
	double fx = std::max(0.0, x0 - x1);
    
	if (!grad.empty()) {
		if (!fixed[ix0]) grad[var_indices[ix0]] = (fx == 0.0? 0 :  1);
        if (!fixed[ix1]) grad[var_indices[ix1]] = (fx == 0.0? 0 : -1);
	}
    
	return  fx;
}


//double constraints::leq(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
//                       std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
//{
//    //- derivative (x1 - x0)^2
//    
//	int ix0 = indices[0];
//	int ix1 = indices[1];
//    
//	double x0 = fixed[ix0] ? fixed_params[fixed_indices[ix0]] : x[var_indices[ix0]];
//	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];
//    
//	double fx = std::min(0.0, x1 - x0);
//	double fx2 = fx *fx;
//    
//	if (!grad.empty()) {
//		if (!fixed[ix0]) grad[var_indices[ix0]] = (fx == 0.0? 0 : 2*(x0 - x1));
//        if (!fixed[ix1]) grad[var_indices[ix1]] = (fx == 0.0? 0 : 2*(x1 - x0));
//	}
//    
//	return  fx2;
//}


/*
 LT: x0 < x1
 */

double constraints::lt(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
                       std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
    //- derivative (x1 - x0)^2
    
	int ix0 = indices[0];
	int ix1 = indices[1];
    
    double tol = 0.1; //- tolerance
    
	double x0 = fixed[ix0] ? fixed_params[fixed_indices[ix0]] : x[var_indices[ix0]];
	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];
    
	double fx = std::max(0.0, x0 + tol - x1);
    
	if (!grad.empty()) {
		if (!fixed[ix0]) grad[var_indices[ix0]] = (fx == 0.0? 0 :  1);
        if (!fixed[ix1]) grad[var_indices[ix1]] = (fx == 0.0? 0 : -1);
	}
    
	return  fx;
}



//double constraints::lt(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
//                        std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
//{
//    //- derivative (x1 - x0)^2
//    
//	int ix0 = indices[0];
//	int ix1 = indices[1];
//    
//    double tol = 0.1; //- tolerance
//    
//	double x0 = fixed[ix0] ? fixed_params[fixed_indices[ix0]] : x[var_indices[ix0]];
//	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];
//    
//	double fx = std::min(0.0, x1 - (x0 + tol));
//	double fx2 = fx *fx;
//    
//	if (!grad.empty()) {
//		if (!fixed[ix0]) grad[var_indices[ix0]] = (fx == 0.0? 0 : 2*(x0 - x1 + tol));
//        if (!fixed[ix1]) grad[var_indices[ix1]] = (fx == 0.0? 0 : 2*(x1 - x0 - tol));
//	}
//    
//	return  fx2;
//}


/*
 DR: circle cirlce
 */
double constraints::dr(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
                       std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int ix1 = indices[0], iy1 = indices[1], ir1 = indices[2];
	int ix2 = indices[3], iy2 = indices[4], ir2 = indices[5];
    
	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];
	double y1 = fixed[iy1] ? fixed_params[fixed_indices[iy1]] : x[var_indices[iy1]];
	double r1 = fixed[ir1] ? fixed_params[fixed_indices[ir1]] : x[var_indices[ir1]];
	double x2 = fixed[ix2] ? fixed_params[fixed_indices[ix2]] : x[var_indices[ix2]];
	double y2 = fixed[iy2] ? fixed_params[fixed_indices[iy2]] : x[var_indices[iy2]];
	double r2 = fixed[ir2] ? fixed_params[fixed_indices[ir2]] : x[var_indices[ir2]];
    
	double dx = (x1 - x2);
	double dy = (y1 - y2);
	double r12 = (r1 + r2);
	double fx = std::min(0.0, dx*dx + dy*dy - r12*r12);
	double fx2 = fx *fx;
    
	if (!grad.empty()) {
		if (!fixed[ix1]) grad[var_indices[ix1]] = (fx == 0.0? 0 : 4 * dx * (fx));
		if (!fixed[iy1]) grad[var_indices[iy1]] = (fx == 0.0? 0 : 4 * dy * (fx));
		if (!fixed[ir1]) grad[var_indices[ir1]] = (fx == 0.0? 0 : -4 * r12 * (fx));
		if (!fixed[ix2]) grad[var_indices[ix2]] = (fx == 0.0? 0 : -4 * dx * (fx));
		if (!fixed[iy2]) grad[var_indices[iy2]] = (fx == 0.0? 0 : -4 * dy * (fx));
		if (!fixed[ir2]) grad[var_indices[ir2]] = (fx == 0.0? 0 : -4 * r12 * (fx));
	}
    
	return  fx2;
}



/*
 DC: circle cirlce
 */
double constraints::dc(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
                       std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int ix1 = indices[0], iy1 = indices[1], ir1 = indices[2];
	int ix2 = indices[3], iy2 = indices[4], ir2 = indices[5];
    
    double tol = 0.1;
    
	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];
	double y1 = fixed[iy1] ? fixed_params[fixed_indices[iy1]] : x[var_indices[iy1]];
	double r1 = fixed[ir1] ? fixed_params[fixed_indices[ir1]] : x[var_indices[ir1]];
	double x2 = fixed[ix2] ? fixed_params[fixed_indices[ix2]] : x[var_indices[ix2]];
	double y2 = fixed[iy2] ? fixed_params[fixed_indices[iy2]] : x[var_indices[iy2]];
	double r2 = fixed[ir2] ? fixed_params[fixed_indices[ir2]] : x[var_indices[ir2]];
    
	double dx = (x1 - x2);
	double dy = (y1 - y2);
	double r12 = (r1 + r2);
	double fx = std::min(0.0, dx*dx + dy*dy - tol - r12*r12);
	double fx2 = fx *fx;
    
	if (!grad.empty()) {
		if (!fixed[ix1]) grad[var_indices[ix1]] = (fx == 0.0? 0 : 4 * dx * (fx));
		if (!fixed[iy1]) grad[var_indices[iy1]] = (fx == 0.0? 0 : 4 * dy * (fx));
		if (!fixed[ir1]) grad[var_indices[ir1]] = (fx == 0.0? 0 : -4 * r12 * (fx));
		if (!fixed[ix2]) grad[var_indices[ix2]] = (fx == 0.0? 0 : -4 * dx * (fx));
		if (!fixed[iy2]) grad[var_indices[iy2]] = (fx == 0.0? 0 : -4 * dy * (fx));
		if (!fixed[ir2]) grad[var_indices[ir2]] = (fx == 0.0? 0 : -4 * r12 * (fx));
	}
    
	return  fx2;
}


/*
 EC: circle cirlce
*/
double constraints::ec(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
	std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices) 
{
	int ix1 = indices[0], iy1 = indices[1], ir1 = indices[2];
	int ix2 = indices[3], iy2 = indices[4], ir2 = indices[5];

	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];
	double y1 = fixed[iy1] ? fixed_params[fixed_indices[iy1]] : x[var_indices[iy1]];
	double r1 = fixed[ir1] ? fixed_params[fixed_indices[ir1]] : x[var_indices[ir1]];
	double x2 = fixed[ix2] ? fixed_params[fixed_indices[ix2]] : x[var_indices[ix2]];
	double y2 = fixed[iy2] ? fixed_params[fixed_indices[iy2]] : x[var_indices[iy2]];
	double r2 = fixed[ir2] ? fixed_params[fixed_indices[ir2]] : x[var_indices[ir2]];

	double dx = (x1 - x2);
	double dy = (y1 - y2);
	double r12 = (r1 + r2);
	double fx = dx*dx + dy*dy - r12*r12;
	double fx2 = fx *fx;

	if (!grad.empty()) {
		if (!fixed[ix1]) grad[var_indices[ix1]] = 4 * dx * (fx);
		if (!fixed[iy1]) grad[var_indices[iy1]] = 4 * dy * (fx);
		if (!fixed[ir1]) grad[var_indices[ir1]] = -4 * r12 * (fx);
		if (!fixed[ix2]) grad[var_indices[ix2]] = -4 * dx * (fx);
		if (!fixed[iy2]) grad[var_indices[iy2]] = -4 * dy * (fx);
		if (!fixed[ir2]) grad[var_indices[ir2]] = -4 * r12 * (fx);
	}

	return  fx2;
}

/*
tangent_internal: circle circle
*/
double constraints::tangent_internal(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
	std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int ix1 = indices[0], iy1 = indices[1], ir1 = indices[2];
	int ix2 = indices[3], iy2 = indices[4], ir2 = indices[5];

	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];
	double y1 = fixed[iy1] ? fixed_params[fixed_indices[iy1]] : x[var_indices[iy1]];
	double r1 = fixed[ir1] ? fixed_params[fixed_indices[ir1]] : x[var_indices[ir1]];
	double x2 = fixed[ix2] ? fixed_params[fixed_indices[ix2]] : x[var_indices[ix2]];
	double y2 = fixed[iy2] ? fixed_params[fixed_indices[iy2]] : x[var_indices[iy2]];
	double r2 = fixed[ir2] ? fixed_params[fixed_indices[ir2]] : x[var_indices[ir2]];

	double dx = (x1 - x2);
	double dy = (y1 - y2);
	double r12 = (r1 - r2);
	double r21 = (r2 - r1);
	double d = sqrt(dx*dx + dy*dy);
	double fx = d - r21;
	double fx2 = fx *fx;

	if (!grad.empty()) {
		if (!fixed[ix1]) grad[var_indices[ix1]] = (2 * dx * (fx)) / d;
		if (!fixed[iy1]) grad[var_indices[iy1]] = (2 * dy * (fx)) / d;
		if (!fixed[ir1]) grad[var_indices[ir1]] = 2 * r21 * (fx);
		if (!fixed[ix2]) grad[var_indices[ix2]] = (-2 * dx * (fx)) / d;
		if (!fixed[iy2]) grad[var_indices[iy2]] = (-2 * dy * (fx)) / d;
		if (!fixed[ir2]) grad[var_indices[ir2]] = -2 * r21 * (fx);
	}

	return  fx2;
}


/*
 non_tangent_internal: circle circle
 */
double constraints::non_tangent_internal(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
                        std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int ix1 = indices[0], iy1 = indices[1], ir1 = indices[2];
	int ix2 = indices[3], iy2 = indices[4], ir2 = indices[5];
    
    double tol = 0.1;
    
	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];
	double y1 = fixed[iy1] ? fixed_params[fixed_indices[iy1]] : x[var_indices[iy1]];
	double r1 = fixed[ir1] ? fixed_params[fixed_indices[ir1]] : x[var_indices[ir1]];
	double x2 = fixed[ix2] ? fixed_params[fixed_indices[ix2]] : x[var_indices[ix2]];
	double y2 = fixed[iy2] ? fixed_params[fixed_indices[iy2]] : x[var_indices[iy2]];
	double r2 = fixed[ir2] ? fixed_params[fixed_indices[ir2]] : x[var_indices[ir2]];
    
	double dx = (x1 - x2);
	double dy = (y1 - y2);
	double r12 = (r1 - r2);
	double r21 = (r2 - r1);
	double d = sqrt(dx*dx + dy*dy);
	double fx = std::max(0.0, d + tol - r21);
	double fx2 = fx *fx;
    
    if(d == 0.0 && fx > 0.0) d = 0.01; //- prevent divide-by-zero error
    
	if (!grad.empty()) {
		if (!fixed[ix1]) grad[var_indices[ix1]] = (fx == 0.0? 0 : (2 * dx * (fx)) / d);
		if (!fixed[iy1]) grad[var_indices[iy1]] = (fx == 0.0? 0 : (2 * dy * (fx)) / d);
		if (!fixed[ir1]) grad[var_indices[ir1]] = (fx == 0.0? 0 : 2 * r21 * (fx));
		if (!fixed[ix2]) grad[var_indices[ix2]] = (fx == 0.0? 0 : (-2 * dx * (fx)) / d);
		if (!fixed[iy2]) grad[var_indices[iy2]] = (fx == 0.0? 0 : (-2 * dy * (fx)) / d);
		if (!fixed[ir2]) grad[var_indices[ir2]] = (fx == 0.0? 0 : -2 * r21 * (fx));
	}
    
	return  fx2;
}



/*
Collinear: point line
*/
double constraints::collinear(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
	std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int ix0 = indices[0], iy0 = indices[1];// point
	int ix1 = indices[2], iy1 = indices[3], ix2 = indices[4], iy2 = indices[5]; // line

	double x0 = fixed[ix0] ? fixed_params[fixed_indices[ix0]] : x[var_indices[ix0]];// point
	double y0 = fixed[iy0] ? fixed_params[fixed_indices[iy0]] : x[var_indices[iy0]];
	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];// line
	double y1 = fixed[iy1] ? fixed_params[fixed_indices[iy1]] : x[var_indices[iy1]];
	double x2 = fixed[ix2] ? fixed_params[fixed_indices[ix2]] : x[var_indices[ix2]];
	double y2 = fixed[iy2] ? fixed_params[fixed_indices[iy2]] : x[var_indices[iy2]];

	double dx = x2 - x1;
	double dy = y2 - y1;
	double d2 = dx*dx + dy*dy;
	double d = sqrt(d2);
	double area = x0*(y1 - y2) + x1*(y2 - y0) + x2*(y0 - y1);
	double fx2 = area * area;

	if (!grad.empty()) {
		if (!fixed[ix0]) grad[var_indices[ix0]] = 2 * (y1 - y2) * area;
		if (!fixed[ix1]) grad[var_indices[ix1]] = 2 * (y2 - y0) * area;
		if (!fixed[ix2]) grad[var_indices[ix2]] = 2 * (y0 - y1) * area;
		if (!fixed[iy0]) grad[var_indices[iy0]] = 2 * (x2 - x1) * area;
		if (!fixed[iy1]) grad[var_indices[iy1]] = 2 * (x0 - x2) * area;
		if (!fixed[iy2]) grad[var_indices[iy2]] = 2 * (x1 - x0) * area;
	}

	return fx2;
}


/*
Perpendicular: line line
*/
double constraints::perpendicular(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
	std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int il1x1 = indices[0], il1y1 = indices[1], il1x2 = indices[2], il1y2 = indices[3]; 
	int il2x1 = indices[4], il2y1 = indices[5], il2x2 = indices[6], il2y2 = indices[7]; 

	double l1p1x = fixed[il1x1] ? fixed_params[fixed_indices[il1x1]] : x[var_indices[il1x1]];
	double l1p2x = fixed[il1x2] ? fixed_params[fixed_indices[il1x2]] : x[var_indices[il1x2]];
	double l2p1x = fixed[il2x1] ? fixed_params[fixed_indices[il2x1]] : x[var_indices[il2x1]];
	double l2p2x = fixed[il2x2] ? fixed_params[fixed_indices[il2x2]] : x[var_indices[il2x2]];

	double l1p1y = fixed[il1y1] ? fixed_params[fixed_indices[il1y1]] : x[var_indices[il1y1]];
	double l1p2y = fixed[il1y2] ? fixed_params[fixed_indices[il1y2]] : x[var_indices[il1y2]];
	double l2p1y = fixed[il2y1] ? fixed_params[fixed_indices[il2y1]] : x[var_indices[il2y1]];
	double l2p2y = fixed[il2y2] ? fixed_params[fixed_indices[il2y2]] : x[var_indices[il2y2]];

	double dx1 = (l1p1x - l1p2x);
	double dy1 = (l1p1y - l1p2y);
	double dx2 = (l2p1x - l2p2x);
	double dy2 = (l2p1y - l2p2y);
	double fx = (dx1*dx2 + dy1*dy2);
	double fx2 = fx*fx;


	if (!grad.empty()) {
		if (!fixed[il1x1]) grad[var_indices[il1x1]] = 2 * dx2 * fx;
		if (!fixed[il1y1]) grad[var_indices[il1y1]] = 2 * dy2 * fx;
		if (!fixed[il1x2]) grad[var_indices[il1x2]] = -2 * dx2 * fx;
		if (!fixed[il1y2]) grad[var_indices[il1y2]] = -2 * dy2 * fx;
		if (!fixed[il2x1]) grad[var_indices[il2x1]] = 2 * dx1 * fx;
		if (!fixed[il2y1]) grad[var_indices[il2y1]] = 2 * dy1 * fx;
		if (!fixed[il2x2]) grad[var_indices[il2x2]] = -2 * dx1 * fx;
		if (!fixed[il2y2]) grad[var_indices[il2y2]] = -2 * dy1 * fx;
	}

	return fx2;
}


/*
Parallel: line line
*/
double constraints::parallel(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
	std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int il1x1 = indices[0], il1y1 = indices[1], il1x2 = indices[2], il1y2 = indices[3];
	int il2x1 = indices[4], il2y1 = indices[5], il2x2 = indices[6], il2y2 = indices[7];

	double l1p1x = fixed[il1x1] ? fixed_params[fixed_indices[il1x1]] : x[var_indices[il1x1]];
	double l1p2x = fixed[il1x2] ? fixed_params[fixed_indices[il1x2]] : x[var_indices[il1x2]];
	double l2p1x = fixed[il2x1] ? fixed_params[fixed_indices[il2x1]] : x[var_indices[il2x1]];
	double l2p2x = fixed[il2x2] ? fixed_params[fixed_indices[il2x2]] : x[var_indices[il2x2]];

	double l1p1y = fixed[il1y1] ? fixed_params[fixed_indices[il1y1]] : x[var_indices[il1y1]];
	double l1p2y = fixed[il1y2] ? fixed_params[fixed_indices[il1y2]] : x[var_indices[il1y2]];
	double l2p1y = fixed[il2y1] ? fixed_params[fixed_indices[il2y1]] : x[var_indices[il2y1]];
	double l2p2y = fixed[il2y2] ? fixed_params[fixed_indices[il2y2]] : x[var_indices[il2y2]];

	double dx1 = (l1p1x - l1p2x);
	double dy1 = (l1p1y - l1p2y);
	double dx2 = (l2p1x - l2p2x);
	double dy2 = (l2p1y - l2p2y);
	double fx = (dx1*dy2 - dy1*dx2);
	double fx2 = fx*fx;

	if (!grad.empty()) {
		if (!fixed[il1x1]) grad[var_indices[il1x1]] = 2 * dy2 * fx;
		if (!fixed[il1y1]) grad[var_indices[il1y1]] = -2 * dx2 * fx;
		if (!fixed[il1x2]) grad[var_indices[il1x2]] = -2 * dy2 * fx;
		if (!fixed[il1y2]) grad[var_indices[il1y2]] = 2 * dx2 * fx;
		if (!fixed[il2x1]) grad[var_indices[il2x1]] = -2 * dy1 * fx;
		if (!fixed[il2y1]) grad[var_indices[il2y1]] = 2 * dx1 * fx;
		if (!fixed[il2x2]) grad[var_indices[il2x2]] = 2 * dy1 * fx;
		if (!fixed[il2y2]) grad[var_indices[il2y2]] = -2 * dx1 * fx;
	}

	return fx2;
}

/*
Coincident: point circle
*/
double constraints::coincident(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
	std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int ix0 = indices[0], iy0 = indices[1];// point
	int ix1 = indices[2], iy1 = indices[3], ir1 = indices[4]; // circle

	double x0 = fixed[ix0] ? fixed_params[fixed_indices[ix0]] : x[var_indices[ix0]]; // point
	double y0 = fixed[iy0] ? fixed_params[fixed_indices[iy0]] : x[var_indices[iy0]];
	double x1 = fixed[ix1] ? fixed_params[fixed_indices[ix1]] : x[var_indices[ix1]];// circle
	double y1 = fixed[iy1] ? fixed_params[fixed_indices[iy1]] : x[var_indices[iy1]];
	double r1 = fixed[ir1] ? fixed_params[fixed_indices[ir1]] : x[var_indices[ir1]];

	double dx = (x0 - x1);
	double dy = (y0 - y1);
	double fx = dx*dx + dy*dy - r1*r1;
	double fx2 = fx*fx;

	if (!grad.empty()) {
		if (!fixed[ix0]) grad[var_indices[ix0]] = 4 * dx * fx;
		if (!fixed[iy0]) grad[var_indices[iy0]] = 4 * dy * fx;
		if (!fixed[ix1]) grad[var_indices[ix1]] = -4 * dx * fx;
		if (!fixed[iy1]) grad[var_indices[iy1]] = -4 * dy * fx;
		if (!fixed[ir1]) grad[var_indices[ir1]] = -4 * r1 * fx;
	}

	return fx2;
}

/*
Angle: line line
*/
double constraints::angle(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
	std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int il1x1 = indices[0], il1y1 = indices[1], il1x2 = indices[2], il1y2 = indices[3];
	int il2x1 = indices[4], il2y1 = indices[5], il2x2 = indices[6], il2y2 = indices[7];
	int iangl = indices[8];

	double l1p1x = fixed[il1x1] ? fixed_params[fixed_indices[il1x1]] : x[var_indices[il1x1]];
	double l1p2x = fixed[il1x2] ? fixed_params[fixed_indices[il1x2]] : x[var_indices[il1x2]];
	double l2p1x = fixed[il2x1] ? fixed_params[fixed_indices[il2x1]] : x[var_indices[il2x1]];
	double l2p2x = fixed[il2x2] ? fixed_params[fixed_indices[il2x2]] : x[var_indices[il2x2]];

	double l1p1y = fixed[il1y1] ? fixed_params[fixed_indices[il1y1]] : x[var_indices[il1y1]];
	double l1p2y = fixed[il1y2] ? fixed_params[fixed_indices[il1y2]] : x[var_indices[il1y2]];
	double l2p1y = fixed[il2y1] ? fixed_params[fixed_indices[il2y1]] : x[var_indices[il2y1]];
	double l2p2y = fixed[il2y2] ? fixed_params[fixed_indices[il2y2]] : x[var_indices[il2y2]];

	double angl = fixed[iangl] ? fixed_params[fixed_indices[iangl]] : x[var_indices[iangl]];

	double dx1 = (l1p2x - l1p1x);
	double dy1 = (l1p2y - l1p1y);
	double dx2 = (l2p2x - l2p1x);
	double dy2 = (l2p2y - l2p1y);
	double a = atan2(dy1, dx1) + angl;
	double ca = cos(a);
	double sa = sin(a);
	double x2 = dx2*ca + dy2*sa;
	double y2 = -dx2*sa + dy2*ca;
	double fx =  atan2(y2, x2);
	double fx2 = fx*fx;

	if (!grad.empty()) {
		double r2 = dx1*dx1 + dy1*dy1;
		if (!fixed[il1x1]) grad[var_indices[il1x1]] = 2 * -dy1 / r2;
		if (!fixed[il1y1]) grad[var_indices[il1y1]] = 2 * dx1 / r2;
		if (!fixed[il1x2]) grad[var_indices[il1x2]] = 2 * dy1 / r2;
		if (!fixed[il1y2]) grad[var_indices[il1y2]] = 2 * -dx1 / r2;

		r2 = dx2*dx2 + dy2*dy2;
		dx2 = -y2 / r2;
		dy2 = x2 / r2;
		if (!fixed[il2x1]) grad[var_indices[il2x1]] = (-2*ca*dx2 + sa*dy2);
		if (!fixed[il2y1]) grad[var_indices[il2y1]] = (-2*sa*dx2 - ca*dy2);
		if (!fixed[il2x2]) grad[var_indices[il2x2]] = (2*ca*dx2 - sa*dy2);
		if (!fixed[il2y2]) grad[var_indices[il2y2]] = (2*sa*dx2 + ca*dy2);
		if (!fixed[iangl]) grad[var_indices[iangl]] = -1;
	}

	return fx2;
}

/*
Diff vc = va - vb 
*/
double constraints::diff(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
	std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int ivc = indices[0], iva = indices[1], ivb = indices[2];
	
	double vc = fixed[ivc] ? fixed_params[fixed_indices[ivc]] : x[var_indices[ivc]];
	double va = fixed[iva] ? fixed_params[fixed_indices[iva]] : x[var_indices[iva]];
	double vb = fixed[ivb] ? fixed_params[fixed_indices[ivb]] : x[var_indices[ivb]];

	double fx = vc - (va - vb);
	double fx2 = fx * fx;

	if (!grad.empty()) {
		if (!fixed[ivc]) grad[var_indices[ivc]] = 2 * fx;
		if (!fixed[iva]) grad[var_indices[iva]] = -2 * fx;
		if (!fixed[ivb]) grad[var_indices[ivb]] = 2 * fx;
	}

	return fx2;
}



/*
f1 =    sqrt((x1-x2)^2 + (y1-y2)^2) * (x4-x3) +  sqrt((x3-x4)^2 + (y3-y4)^2) * (y2-y1)

f2 =    sqrt((x1-x2)^2 + (y1-y2)^2) * (y4-y3) -  sqrt((x3-x4)^2 + (y3-y4)^2) * (x2-x1)


f = f1^2 + f2^2
*/
double constraints::angle90cc(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
	std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int il1x1 = indices[0], il1y1 = indices[1], il1x2 = indices[2], il1y2 = indices[3];
	int il2x1 = indices[4], il2y1 = indices[5], il2x2 = indices[6], il2y2 = indices[7];

	double x1 = fixed[il1x1] ? fixed_params[fixed_indices[il1x1]] : x[var_indices[il1x1]];
	double x2 = fixed[il1x2] ? fixed_params[fixed_indices[il1x2]] : x[var_indices[il1x2]];
	double x3 = fixed[il2x1] ? fixed_params[fixed_indices[il2x1]] : x[var_indices[il2x1]];
	double x4 = fixed[il2x2] ? fixed_params[fixed_indices[il2x2]] : x[var_indices[il2x2]];

	double y1 = fixed[il1y1] ? fixed_params[fixed_indices[il1y1]] : x[var_indices[il1y1]];
	double y2 = fixed[il1y2] ? fixed_params[fixed_indices[il1y2]] : x[var_indices[il1y2]];
	double y3 = fixed[il2y1] ? fixed_params[fixed_indices[il2y1]] : x[var_indices[il2y1]];
	double y4 = fixed[il2y2] ? fixed_params[fixed_indices[il2y2]] : x[var_indices[il2y2]];

	double dx1 = (x1 - x2);
	double dx2 = (x3 - x4);
	double dy1 = (y1 - y2);
	double dy2 = (y3 - y4);
	double d12 = sqrt(dx1*dx1 + dy1*dy1);
	double d34 = sqrt(dx2*dx2 + dy2*dy2);

	double f1 = (-dx2 * d12) + (-dy1 * d34 );
	double f2 = (-dy2 * d12) - (-dx1 * d34);
	double fx = f1*f1 + f2*f2;

	if (!grad.empty()) {
		if (!fixed[il1x1]) grad[var_indices[il1x1]] = (-2*dx2 * dx1 * (-dx2*d12 -dy1*d34)) /d12 + 2 * ((-dy2*dx1 / d12) + d34)*(-dy2*d12 + dx1*d34);
		if (!fixed[il1y1]) grad[var_indices[il1y1]] = 2 * (((-dx2*dy1) / d12) - d34)*((-dx2)*d12 + (-dy1)*d34) + ((2 * dy1*(-dy2)*(dx1*d34 -dy2*d12)) / d12);
		if (!fixed[il1x2]) grad[var_indices[il1x2]] = 2 * (-((dx1*(-dy2)) / d12) - d34)*(-dy2*d12 + dx1 *d34) - ((2 * dx1*(-dx2) * (-dx2*d12 -dy1 * d34)) / (d12));
		if (!fixed[il1y2]) grad[var_indices[il1y2]] = 2 * ((dx2*dy1) / d12 + d34)*((-dx2)*d12 -dy1 * d34) - (2 * dy1*(-dy2)*(dx1*d34 - dy2 *d12)) / d12;
		if (!fixed[il2x1]) grad[var_indices[il2x1]] = 2 * (((dx2*(-dy1)) / d34) - d12)*(-dx2*d12 -dy1 * d34) - ((2 * (-dx1)*dx2 * (dx1*d34 - dy2 * d12)) / (d34));
		if (!fixed[il2y1]) grad[var_indices[il2y1]] = (2 * (-dy1)*dy2*((-dx2)*d12 -dy1 * d34)) / d34 + 2*((dx1*dy2) / d34 - d12)*(dx1*d34 -dy2 * d12);
		if (!fixed[il2x2]) grad[var_indices[il2x2]] = 2 * (d12+ ((dx2*dy1) / d34) )*(-dx2*d12 -dy1 * d34) + (2 * (-dx1)*dx2 * (dx1*d34 -dy2 * d12)) / (d34);
		if (!fixed[il2y2]) grad[var_indices[il2y2]] = 2 * ((-dx1*dy2) / d34 + d12)*(dx1*d34 -dy2 *d12) - (2 * (-dy1)*(dy2)*((-dx2)*d12 -dy1 * d34)) / d34;
	}

	return fx;
}


/*
f1 =    sqrt((x1-x2)^2 + (y1-y2)^2) * (x4-x3)
-  sqrt((x3-x4)^2 + (y3-y4)^2) * (x2-x1)

f2 =    sqrt((x1-x2)^2 + (y1-y2)^2) * (y4-y3)
-  sqrt((x3-x4)^2 + (y3-y4)^2) * (y2-y1)
*/
double constraints::angle_eq(const std::vector<double> &x, std::vector<double> &grad, std::vector<int>& indices,
	std::vector<bool>& fixed, std::vector<double>& fixed_params, std::map<int, int>& fixed_indices, std::map<int, int>& var_indices)
{
	int il1x1 = indices[0], il1y1 = indices[1], il1x2 = indices[2], il1y2 = indices[3];
	int il2x1 = indices[4], il2y1 = indices[5], il2x2 = indices[6], il2y2 = indices[7];

	double x1 = fixed[il1x1] ? fixed_params[fixed_indices[il1x1]] : x[var_indices[il1x1]];
	double x2 = fixed[il1x2] ? fixed_params[fixed_indices[il1x2]] : x[var_indices[il1x2]];
	double x3 = fixed[il2x1] ? fixed_params[fixed_indices[il2x1]] : x[var_indices[il2x1]];
	double x4 = fixed[il2x2] ? fixed_params[fixed_indices[il2x2]] : x[var_indices[il2x2]];

	double y1 = fixed[il1y1] ? fixed_params[fixed_indices[il1y1]] : x[var_indices[il1y1]];
	double y2 = fixed[il1y2] ? fixed_params[fixed_indices[il1y2]] : x[var_indices[il1y2]];
	double y3 = fixed[il2y1] ? fixed_params[fixed_indices[il2y1]] : x[var_indices[il2y1]];
	double y4 = fixed[il2y2] ? fixed_params[fixed_indices[il2y2]] : x[var_indices[il2y2]];

	double dx1 = (x1 - x2);
	double dx2 = (x3 - x4);
	double dy1 = (y1 - y2);
	double dy2 = (y3 - y4);
	double d12 = sqrt(dx1*dx1 + dy1*dy1);
	double d34 = sqrt(dx2*dx2 + dy2*dy2);

	double f1 = (-dx2 * d12) - (-dx1 * d34);
	double f2 = (-dy2 * d12) - (-dy1 * d34 );
	double fx = f1*f1 + f2*f2;

	if (!grad.empty()) {
		if (!fixed[il1x1]) grad[var_indices[il1x1]] = 2 * (dx1*(-dx2) / d12 + d34)*(-dx2*d12 + dx1*d34) + (2 * dx1*(-dy2) * (-dy2*d12 + dy1*d34)) / d12;
		if (!fixed[il1y1]) grad[var_indices[il1y1]] = (2 * (-dx2)*dy1 * (-dx2*d12 + dx1*d34)) / d12 + 2*(-dy2 + d12 + dy1*d34)*(-dy2*dy1 / d12 + d34);
		if (!fixed[il1x2]) grad[var_indices[il1x2]] = 2 * (-dx1*(-dx2)/d12 - d34)*(-dx2*d12 + dx1*d34) - (2 * dx1*(-dy2) * (-dy2*d12 + dy1*d34)) / d12;
		if (!fixed[il1y2]) grad[var_indices[il1y2]] = 2 * (-dy2*d12 + dy1*d34)*(-(dy1*(-dy2))/d12 - d34) - (2 * (-dx2)*dy1*(-dx2*d12 + dx1*d34)) / d12;
		if (!fixed[il2x1]) grad[var_indices[il2x1]] = 2 * (-(-dx1)*dx2/d34 - d12)*(-dx2*d12 + dx1*d34) -(2 * dx2*(-dy1) * (-dy2*d12 + dy1*d34)) / d34;
		if (!fixed[il2y1]) grad[var_indices[il2y1]] = 2 * (dy1*dy2/d34 - d12)*(-dy2*d12 + dy1*d34) -(2 * (-dx1)*dy2*(-dx2*d12 + dx1*d34)) / d34;
		if (!fixed[il2x2]) grad[var_indices[il2x2]] = 2 * ((-dx1)*dx2 / d34 + d12)*(-dx2*d12 + dx1*d34) + (2 * dx2*(-dy1) * (-dy2*d12 + dy1*d34)) / d34;
		if (!fixed[il2y2]) grad[var_indices[il2y2]] = (2 * (-dx1)*dy2*(-dx2*d12 + dx1*d34))/d34 + 2*(d12 + (-dy1)*dy2/d34)*(-dy2*d12 + dy1*d34);
	}

	return fx;
}
