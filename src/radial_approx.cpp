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


#include "dmp/radial_approx.h"
#include<stdio.h>
using namespace Eigen;
using namespace std;

namespace dmp{

RadialApprox::RadialApprox(int num_bases, double base_width, double alpha)
{
	n_bases = num_bases;
	features = new double[n_bases];
	centers = new double[n_bases];
	widths = new double[n_bases];
	weights.resize(n_bases);
	for(int i=0; i<n_bases; i++){
		features[i] = 0;
		centers[i] = exp((-alpha*i)/n_bases);
		widths[i] = base_width * (1/exp((-alpha*i)/n_bases));
	}
}


RadialApprox::RadialApprox(const vector<double> &w, double base_width, double alpha)
{
	weights = w;
	n_bases = w.size();
	features = new double[n_bases];
		centers = new double[n_bases];
		widths = new double[n_bases];
		for(int i=0; i<n_bases; i++){
			features[i] = 0;
			centers[i] = ((double)i)/((double)n_bases);  //exp((-alpha*i)/n_bases);
			widths[i] = base_width; //base_width * exp((-alpha*i)/n_bases);
		}
}


RadialApprox::~RadialApprox()
{
	delete[] features;
	delete[] centers;
	delete[] widths;
}


double RadialApprox::evalAt(double x)
{
	calcFeatures(x);

	double wsum = 0;
	for(int i=0; i<n_bases; i++){
		wsum += features[i] * weights[i];
	}
	return wsum;
}


void RadialApprox::leastSquaresWeights(double *X, double *Y, int n_pts)
{
	MatrixXd D_mat = MatrixXd(n_pts,n_bases);
	MatrixXd Y_mat = MatrixXd(n_pts,1);

	//Calculate the design matrix
	for(int i=0; i<n_pts; i++){
		Y_mat(i,0) = Y[i];
		calcFeatures(X[i]);
		for(int j=0; j<n_bases; j++){
			D_mat(i,j) = features[j];
		}
	}

	//Calculate the least squares weights via projection onto the basis functions
	MatrixXd w = pseudoinverse(D_mat.transpose() * D_mat) * D_mat.transpose() * Y_mat;
	for(int i=0; i<n_bases; i++){
		weights[i] = w(i,0);
	}
}


void RadialApprox::calcFeatures(double x)
{
	double sum = 0;
	for(int i=0; i<n_bases; i++){
		features[i] = exp(-widths[i]*(x-centers[i])*(x-centers[i]));
		sum += features[i];
	}
	for(int i=0; i<n_bases; i++){
		features[i] /= sum;
	}
}


MatrixXd RadialApprox::pseudoinverse(MatrixXd mat){
	//Numpy uses 1e-15 by default.  I use 1e-10 just to be safe.
	double precisionCutoff = 1e-10;

	//Compute the SVD of the matrix
	JacobiSVD<MatrixXd> svd(mat, ComputeThinU | ComputeThinV);
	MatrixXd U = svd.matrixU();
	MatrixXd V = svd.matrixV();
	MatrixXd S = svd.singularValues();

	//Psuedoinvert the diagonal matrix of singular values
	MatrixXd S_plus = MatrixXd::Zero(n_bases, n_bases);
	for(int i=0; i<n_bases; i++){
		if(S(i) > precisionCutoff){  //Cutoff to avoid huge inverted values for numerical stability
			S_plus(i,i) = 1.0/S(i);
		}
	}

	//Compute psuedoinverse of orginal matrix
	return V * S_plus * U.transpose();
}

}


