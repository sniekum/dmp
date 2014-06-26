/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2012, Scott Niekum
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Robert Bosch nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************/

/**
  * \author Scott Niekum
  */


#include "dmp/linear_approx.h"
#include<stdio.h>
using namespace std;

namespace dmp{
    
    
bool sort_pt_pair(const pt_pair& left, const pt_pair& right)
{
    return left.first < right.first;
}    


LinearApprox::LinearApprox()
{
	n_bases = 0;
}


LinearApprox::LinearApprox(std::vector<double> X, std::vector<double> Y)
{
        n_bases = X.size();
        
        // Read in points and sort them by X-value
        for(int i=0; i<n_bases; i++){
            points.push_back(pt_pair(X[i], Y[i]));
        }
        
        std::sort(points.begin(), points.end(), sort_pt_pair);
        
}


LinearApprox::~LinearApprox(){};


//Assumes that x is between 0 and 1 inclusive and that fxn is zero at x=0 and x=1 if not specified
//Perform linear interpolation between saved points
double LinearApprox::evalAt(double x)
{
        // If out of bounds, or x=0, return 0
        if(x <= 0.0 || x > 1.0) 
            return 0.0;
    
        //If not points to interp, return 0
        if(n_bases <= 0)
            return 0.0;
                
        // If less than the smallest entry, interp with x=0
        if(x < points[0].first){
            double slope = points[0].second / points[0].first;
            return slope * x;
        }    
        // If greater than largest entry, interp with fxn=0 at x=1
        else if(x > points[n_bases-1].first){    
            double inv_slope = points[n_bases].second / points[n_bases].first; 
            return inv_slope * (1.0 - x);
        }
        // Otherwise, normal interp
        else{
            double curr = 0.0;
            int i = 0;
        
            while(x > curr && i < n_bases){
                i++;
                curr = points[i].first;
            }
            
            double diffx = points[i].first - points[i-1].first;
            double diffy = points[i].second - points[i-1].second;
            double slope = diffy/diffx;
                
            double y_start =  points[i-1].second;
            double x_dist = x - points[i-1].first;
            return y_start + (x_dist * slope);
            
        }
            
            
}

// Nothing to do here, weights are not used
void LinearApprox::leastSquaresWeights(double *X, double *Y, int n_pts){};

}


