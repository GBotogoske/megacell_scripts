#define n_peaks 8
#define draw 8

#include <iostream> 


Double_t gaus_error(Double_t A, Double_t sigma, Double_t A_e, Double_t sigma_e)
{
	Double_t error= sqrt(2*M_PI)*A*sigma*(sqrt(pow(sigma_e/sigma,2)+pow(A_e/A,2)));
	return error;
}

int peak_actual=1;

//função que calcula fatorial
int factorial(int n)
{
	if(n>1)
	{
		return n*factorial(n-1);
	}
	else
		return 1;
}

Double_t gaus_plus(Double_t *x,Double_t *par) //A0, mean0, std0, mean_spe, A1, std1, A2, std2, A3, std3, ... 
{
	Double_t arg = 0;
	Double_t fitval = 0;
	for(int i=0;i<peak_actual;i++)
	{
		if(i==0)
		{
			if (par[2]!=0)
			{
				arg = (x[0] - par[1])/par[2];
			}
			fitval += par[0]*TMath::Exp(-0.5*arg*arg);
		}
		else
		{
			if (par[4]!=0)
			{
				arg = (x[0] - (i*par[3]+par[1]))/(par[3+2*i]);
			}
			fitval += par[3+2*i-1]*TMath::Exp(-0.5*arg*arg);
		}
	}
	return fitval;
}


void fit_hist_charge_3(string filename)
{
	
	TFile *f1 = new TFile(filename.c_str(),"UPDATE"); //abre o .root, atualizando-o caso necessario
	
	TTree *t2 = (TTree*)f1->Get("TChan"); //pega a ttree
	
	//channel channel_variable; //cria uma variavel da estrutura channel
	
	double charge;
	double noise;
	double baseline;
	long long index;
	t2->SetBranchAddress("roiq",&charge);//aponta a variavel criada, para a localização correta na ttree
	t2->SetBranchAddress("filt_noise",&noise);//aponta a variavel criada, para a localização correta na ttree
	t2->SetBranchAddress("filt_baseline",&baseline);//aponta a variavel criada, para a localização correta na ttree
	t2->SetBranchAddress("index",&index);//aponta a variavel criada, para a localização correta na ttree
	
	Int_t entries=(Int_t)t2->GetEntries(); //pega o numero de eventos
	cout<<entries<<endl;
	
	double max_bin=8e3;
	double min_bin=-2e3;
	int n_bin=200;
	
	double size_bin=(max_bin-min_bin)/n_bin;
	
	TH1F *hQ = new TH1F("histCharge","Charge;Charge[ADC];Events",200,-2e3,8e3);
	
	
	//seta ruido maximo permitido
	double noise_max=1.7;
	double baseline_mean=3550;
	double baseline_variation=5;
	
	
	int final_entries=0;
	//preenche o histograma
	for(int i=0;i<entries;i++)
	{
		t2->GetEntry(i);
		if(index==3)
		{
			if(noise<=noise_max)
			{
				//if(abs(baseline-baseline_mean)<baseline_variation)
				//{
					hQ->Fill(charge);
					final_entries++;
				//}
			}
		}
	}
	
	//cria janela para colocar os graficos
	Double_t w = 600;
	Double_t h = 600;
	auto c = new TCanvas("c", "c", w, h);
	c->SetWindowSize(w + (w - c->GetWw()), h + (h - c->GetWh()));
	
	hQ->Draw("");
	
	
	//cria um vetor para definir os parametros do fit:  A0, mean0, std0, mean_spe, std_spe, A1, A2, A3, ...
	Double_t param[2*n_peaks+2];

	float min_fit;
	float max_fit;
	float delta_fit;
	vector<TF1*> gaus_fit(n_peaks); //fit
	
	
	//----------------------------------
	for(;peak_actual<=n_peaks;peak_actual++)
	{
		if(peak_actual==1)
		{
			//baseline --chutes iniciais do ruido
			param[0]=100.0; //A
			param[1]=0.0; //MEAN
			param[2]=100.0; //STD
			
			//região do fit inicial
			delta_fit=200;
			min_fit=-delta_fit;
			max_fit=+delta_fit;	
		}
		else if(peak_actual==2)
		{
			//chute inicial da media e std do spe
			param[3]=700; //mean_spe
			param[5]=150; //std_spe
			
			//todos os picos
			min_fit=-1000;
			max_fit=param[1]+(peak_actual-1)*param[3]+param[5]*sqrt(peak_actual-1);
		}
		else
		{
				//todos os picos
			min_fit=-1000;
			max_fit=param[1]+(peak_actual-1)*param[3]+param[5]*sqrt(peak_actual-1);
			cout<<endl<<max_fit<<endl;
		}
		
		
		//fit total
		if(peak_actual!=1)
			gaus_fit[peak_actual-1] = new TF1("gaus_plus",gaus_plus,min_fit,max_fit,2*peak_actual+2);
		else
			gaus_fit[peak_actual-1] = new TF1("gaus_plus",gaus_plus,min_fit,max_fit,3);	
		
		//coloca nome nas variaveis: A0, mean0, std0, mean_spe, std_spe, A1, A2, A3, ...
		gaus_fit[peak_actual-1]->SetParName(0,"A0");
		gaus_fit[peak_actual-1]->SetParName(1,"Mean0");
		gaus_fit[peak_actual-1]->SetParName(2,"Std0");
	
		gaus_fit[peak_actual-1]->SetParLimits(0,1e2,3e3);	
		gaus_fit[peak_actual-1]->SetParLimits(1,10,300);
		gaus_fit[peak_actual-1]->SetParLimits(2,50,500);
			
		if(peak_actual>=2)
		{
			gaus_fit[peak_actual-1]->SetParName(3,"Mean_spe");
			//gaus_fit[peak_actual-1]->SetParName(4,"Std_spe");
			
			for(int i=1;i<peak_actual;i++)
			{
				gaus_fit[peak_actual-1]->SetParName(3+2*i,Form("std%i",i));
				gaus_fit[peak_actual-1]->SetParName(3+2*i-1,Form("A%i",i));
			}
			
			gaus_fit[peak_actual-1]->SetParLimits(3,100,1000);
			//gaus_fit[peak_actual-1]->SetParLimits(4,50,500);
		//cout<<lambda<<endl;
		
		
			for(int i=1;i<peak_actual;i++)
			{
				gaus_fit[peak_actual-1]->SetParLimits(3+2*i-1,0,3e3);
				gaus_fit[peak_actual-1]->SetParLimits(3+2*i,100,500);
			}
		}
		
		
		gaus_fit[peak_actual-1]->SetParameters(param);

		
		for(int i=0;i<1;i++)
			hQ->Fit(Form("gaus_plus"),"0R");
		
		gaus_fit[peak_actual-1]->GetParameters(param);
			
		if(peak_actual==draw)
		{
			gaus_fit[peak_actual-1]->SetLineColor(1);
			gaus_fit[peak_actual-1]->Draw("SAME");
		
			//colocar a parte de printar os picos aqui
			//....
			//plota cada pico
			vector<TF1*> gaus_peaks(peak_actual);
			
			for(int i=0;i<peak_actual;i++)
			{
				gaus_peaks[i]= new TF1(Form("gaus_peak%i",i),"gaus",-2e3,10e3);
				if(i==0)
				{
					gaus_peaks[i]->SetParameter(0,param[0]);
					gaus_peaks[i]->SetParameter(1,param[1]);
					gaus_peaks[i]->SetParameter(2,param[2]);
				}
				else
				{
					gaus_peaks[i]->SetParameter(0,param[3+2*i-1]);
					gaus_peaks[i]->SetParameter(1,i*param[3]+param[1]);
					gaus_peaks[i]->SetParameter(2,param[3+2*i]);
				}
				gaus_peaks[i]->SetLineColor(i+2);
				gaus_peaks[i]->Draw("SAME");
				
			}
			
			//wait for button
		}
		/*
		if(peak_actual!=n_peaks)
		{
					
			for(int i=0;i<peak_actual;i++)
			{
				delete gaus_peaks[i];
			}
			
			delete gaus_fit;
		}*/
		//delete gaus_fit;
	}
		
	
	//pega os erros calculados
	const Double_t * err;
	err=gaus_fit[n_peaks-1]->GetParErrors();	
	
	//pega o valor do chi2
	Double_t chi2   = gaus_fit[n_peaks-1]->GetChisquare();
	
	//calcula o ganho
	//double gain= (param[3]-param[1])*1e7/1.6;
	
	
	//cria arquivo
	string new_name = "fit.txt";
	ofstream MyFile(new_name.c_str());
	
	MyFile << "A0: " << param[0] << " +- " << err[0] << endl;
	MyFile << "mean0: " << param[1]  << " +- " << err[1] << endl;
	MyFile << "std0: " << param[2] << " +- " << err[2] <<endl;
	MyFile << "mean_spe: " << param[3]  << " +- " << err[3] <<endl;
	
	
	for(int i=1;i<n_peaks;i++)
	{
		MyFile << Form("A%i: ",i) << param[3+2*i-1]  << " +- " << err[3+2*i-1] <<endl;
		MyFile << Form("std%i: ",i) << param[3+2*i]  << " +- " << err[3+2*i] <<endl;
	}
	
	
	cout<<"Chi2: " << chi2 << endl;
		
	MyFile << "Chi2: " << chi2 << endl;

	MyFile.close();

	string new_name_2 = "data.txt";
	ofstream Myfile_2(new_name_2.c_str());
	Myfile_2 << param[0] * sqrt(2*M_PI) * param[2]/(final_entries*size_bin) <<" "<< gaus_error(param[0],param[2],err[0],err[2])/(final_entries*size_bin) << endl;
	
	for(int i=1;i<n_peaks;i++)
	{
		Myfile_2 <<  param[3+2*i-1] * sqrt(2*M_PI) * param[3+2*i]/(final_entries*size_bin) << " "<<gaus_error(param[3+2*i-1],param[3+2*i],err[3+2*i-1],err[3+2*i])/(final_entries*size_bin)<<endl;
	}
	
	Myfile_2.close();
	
	
	
	
	
	
		
	f1->WriteObject(c,"fit","TObject::kOverwrite");
	
}

