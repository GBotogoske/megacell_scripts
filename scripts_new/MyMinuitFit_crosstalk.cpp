#include "TMinuit.h"
#include <TH1F.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <vector>

#define N 5

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
        	Double_t diff = par[2]*y_e[k] - MyFitFunction(k, par);
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
			//printf("%f\n",y_e[i]);
        	x[i]=i;
        }
    	}
        infile.close();
		

	//------------------------------------------------------------------------------
	// Create an instance of the TMinuit class
	TMinuit *minuit = new TMinuit(3);

	// Set the function to be minimized
	minuit->SetFCN(MyChi2Function);

	// Set the starting values and step sizes for the fit parameters

	double lambda=-log(y_e[0]);
	double p=1-y_e[1]/(lambda*y_e[0]);
	double norm=1.2;
	
	Int_t ierr=1;

	Double_t vstart[3] = {p,lambda,norm};
	Double_t step[3] = {1e-6,1e-4,1e-6};
	minuit->mnparm(0, "p", vstart[0], step[0], 0.19, 0.22, ierr);
	minuit->mnparm(1, "lambda", vstart[1], step[1], 0.8, 1, ierr);
	minuit->mnparm(2, "norm", vstart[2], step[2], 0.9, 1.2, ierr);
	

	// Perform the minimization
	minuit->Migrad();

	// Retrieve the results of the fit
	Double_t par[3], err[3];
	minuit->GetParameter(0, par[0], err[0]);
	minuit->GetParameter(1, par[1], err[1]);
	minuit->GetParameter(2, par[2], err[2]);

	
	// Print the results of the fit
	printf("Fit results:\n");
	printf("p = %f +/- %f\n", par[0], err[0]);
	printf("lambda = %f +/- %f\n", par[1], err[1]);
	printf("normalização = %f +/- %f\n", par[2], err[2]);
	

	TCanvas *canvas = new TCanvas("canvas", "My Graph", 1800, 1000);
	
	double y;
	for (Int_t k = 0; k < N; k++)
	{
		y=y_e[k]*par[2];
		//y_e[k]=y_e[k]*par[2];
		y_error[k]=y*sqrt(pow(y_error[k]/y_e[k],2)+pow(err[2]/par[2],2));
		y_e[k]=y_e[k]*par[2];
	}	
	TGraphErrors *graph = new TGraphErrors(N,&x[0],&y_e[0],0, &y_error[0]);
	graph->SetMarkerSize(4);

	// Create a canvas and draw the graph

	graph->SetFillColor(40);
	graph->Draw("AB");
		
	 for (Int_t k = 0; k < N; k++)
	{
		y_f[k]= MyFitFunction(k,par);
	}	

	//vamos estimar os erros do fit
	double erro_1_mais[N], erro_1_menos[N];
	double erro_2_mais[N], erro_2_menos[N];
	Double_t par_0[2], par_1[2],par_2[2], par_3[2];
	
	par_0[0]=par[0];
	par_0[1]=par[1]+err[1];
	
	par_1[0]=par[0];
	par_1[1]=par[1]-err[1];
	
	par_2[0]=par[0]+err[0];
	par_2[1]=par[1];
	
	par_3[0]=par[0]-err[0];
	par_3[1]=par[1];
	
	for (Int_t k = 0; k < N; k++)
	{
		erro_1_mais[k]= MyFitFunction(k,par_0);
		erro_1_menos[k]= MyFitFunction(k,par_1);
		erro_2_mais[k]= MyFitFunction(k,par_2);
		erro_2_menos[k]= MyFitFunction(k,par_3);
		
		erro_1_mais[k]=abs(erro_1_mais[k]-y_f[k]);
		erro_1_menos[k]=abs(erro_1_menos[k]-y_f[k]);
		if(erro_1_menos[k]>erro_1_mais[k])
			erro_1_mais[k]=erro_1_menos[k];
			
		erro_2_mais[k]=abs(erro_2_mais[k]-y_f[k]);
		erro_2_menos[k]=abs(erro_2_menos[k]-y_f[k]);
		if(erro_2_menos[k]>erro_2_mais[k])
			erro_2_mais[k]=erro_2_menos[k];
			
		erro_1_mais[k]=sqrt(erro_1_mais[k]*erro_1_mais[k]+erro_2_mais[k]*erro_2_mais[k]);	
	}	


	TGraphErrors *graph2 = new TGraphErrors(N,&x[0],&y_f[0],0, &erro_1_mais[0]);
	graph2->SetLineColor(2);
	graph->SetFillColor(41);
	graph2->Draw("SAME");
	// Set axis titles
	graph->GetXaxis()->SetTitle("X");
	graph->GetYaxis()->SetTitle("Y");

	auto legend = new TLegend(0.48,0.7,0.98,0.9);
   	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legend->AddEntry(graph,"Experimental data","f");
   	legend->AddEntry(graph2,Form("fitted data -- p=%.2f  -- lambda=%.2f -- kdup=%.2f",par[0],par[1],par[0]/(1-par[0])),"f");
	legend->Draw();
	    graph->GetXaxis()->SetRangeUser(-1, 8);
	// Update the canvas
	canvas->Update();
	canvas->SaveAs(Form("crosstalk.pdf"));
	
	cout<<endl<<chi_2<<endl;
	
  
}


