#include "TMinuit.h"
#include <TH1F.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <vector>

#define N 8

double y_e[N];
double y_error[N];
double y_f[N];
double chi_2;

int factorial(int n) 
{
	int result = 1;
	for(int i=1; i<=n; i++)
	{
	result *= i;
	}
	return result;
}

// Define the function to fit and the chi-squared function
double B(int i, int k)
{
	double val;
	if(i==0 && k==0)
	{
		val=1.0;
	}
	else if(i==0 && k>0)
	{
		val=0.0;
	}
	else
	{
		val=(double)factorial(k-1)/(factorial(i)*factorial(i-1)*factorial(k-i));
	}
	return val;
}

double MyFitFunction(int k, Double_t *par) 
{
	double P=0;
	double p=par[0];
	double lambda = par[1];
	for(int i=0;i<=k;i++)
	{
		P+=B(i,k)*pow(lambda*(1-p),i)*pow(p,k-i);
	}	

	return P*exp(-lambda);
}

void MyChi2Function(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) 
{
	Double_t chi2 = 0.0;
    	for (Int_t k = 0; k < N; k++) 
    	{
        	Double_t diff = y_e[k] - MyFitFunction(k, par);
       	chi2 += diff * diff/y_e[k];
   	}
	f = chi2;
	chi_2=chi2;	
}

void MyMinuitFit_crosstalk() 
{
	
	ifstream infile("data.txt");
 	double x[N];
    	if (infile.is_open()) 
    	{
        double a, b;
        
        for(int i=0;i<N;i++)
        {
        	infile >> y_e[i] >> y_error[i];
        	x[i]=i;
        }
    	}
        infile.close();
		

	//------------------------------------------------------------------------------
	// Create an instance of the TMinuit class
	TMinuit *minuit = new TMinuit(2);

	// Set the function to be minimized
	minuit->SetFCN(MyChi2Function);

	// Set the starting values and step sizes for the fit parameters

	double p=0.27;
	double lambda=1.86;
	
	Int_t ierr=1;

	Double_t vstart[2] = {p,lambda};
	Double_t step[2] = {1e-8,1e-8};
	minuit->mnparm(0, "par0", vstart[0], step[0], 0.26, 0.29, ierr);
	minuit->mnparm(1, "par1", vstart[1], step[1], 1.8, 1.9, ierr);
	

	// Perform the minimization
	minuit->Migrad();

	// Retrieve the results of the fit
	Double_t par[2], err[2];
	minuit->GetParameter(0, par[0], err[0]);
	minuit->GetParameter(1, par[1], err[1]);

	
	// Print the results of the fit
	printf("Fit results:\n");
	printf("par0 = %f +/- %f\n", par[0], err[0]);
	printf("par1 = %f +/- %f\n", par[1], err[1]);

	TCanvas *canvas = new TCanvas("canvas", "My Graph", 800, 600);
	//TGraph *graph = new TGraph();
	TGraphErrors *graph = new TGraphErrors(N,&x[0],&y_e[0],0, &y_error[0]);
	graph->SetMarkerSize(4);

	// Create a canvas and draw the graph

	graph->Draw("AP");
		
	 for (Int_t k = 0; k < N; k++)
	{
		y_f[k]= MyFitFunction(k,par);
	}	

	TGraphErrors *graph2 = new TGraphErrors(N,&x[0],&y_f[0],0, 0);
	graph2->SetLineColor(2);
	graph2->Draw("SAME");
	// Set axis titles
	graph->GetXaxis()->SetTitle("X");
	graph->GetYaxis()->SetTitle("Y");

	auto legend = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legend->AddEntry(graph,"Experimental data","f");
   	legend->AddEntry(graph2,"fitted data","f");
	legend->Draw();
	    graph->GetXaxis()->SetRangeUser(-1, 8);
	// Update the canvas
	canvas->Update();
	
	cout<<endl<<chi_2<<endl;
	
  
}


