void fit_alpha_charge(string filename,int n_channel=0)
{
	TFile *f1 = new TFile(filename.c_str(),"UPDATE"); //abre o .root, atualizando-o caso necessario
	
	TTree *t2 = (TTree*)f1->Get("TChan"); //pega a ttree
	
	//channel channel_variable; //cria uma variavel da estrutura channel
	
	double charge;
	double noise;
	double baseline;
	double fprompt;
	long long index;
	t2->SetBranchAddress("filt_roiq",&charge);//aponta a variavel criada, para a localização correta na ttree
	t2->SetBranchAddress("filt_noise",&noise);//aponta a variavel criada, para a localização correta na ttree
	t2->SetBranchAddress("filt_baseline",&baseline);//aponta a variavel criada, para a localização correta na ttree
	t2->SetBranchAddress("index",&index);//aponta a variavel criada, para a localização correta na ttree
	t2->SetBranchAddress("fprompt",&fprompt);//aponta a variavel criada, para a localização correta na ttree
	
	Int_t entries=(Int_t)t2->GetEntries(); //pega o numero de eventos
	cout<<entries<<endl;
	
	
	//dados para o histograma de carga
	double max_bin=2e6;
	double min_bin=0;
	int n_bin=500;
	
	TH1F *hQ = new TH1F("histCharge","Charge;Charge[ADC*4ns];Events",n_bin,min_bin,max_bin);
	
	double x_min=0;
	double x_max=10e6;
	int n_bin2=1000;
	TH2F *hf = new TH2F("fprompt_cut","fprompt_cut",n_bin2,x_min,x_max,n_bin2,0,1);
	//--------------------------
	
	//varre o arquivo .root
	double noise_max=10;
	double baseline_var=10;
	for(int i=0;i<entries;i++)
	{
		t2->GetEntry(i);
		if(index==n_channel)
		{
			if(noise<=noise_max)
			{
				if(abs(baseline_var-baseline)>0)
				{
					hQ->Fill(charge);
					hf->Fill(charge,fprompt);
				}
			}
		}
	}
	
	
	//janela grafica
	Double_t w = 1800;
	Double_t h = 1000;
	auto c = new TCanvas("c", "c", w, h);
	c->SetWindowSize(w + (w - c->GetWw()), h + (h - c->GetWh()));
	
	double f_cut=0.6;
	double f_y[2]={f_cut,f_cut};
	double f_x[2]={x_min,x_max};
	auto plot1 = new TGraph (2,f_x,f_y );
	hf->SetTitle("my plot;Charge;fprompt");
	hf->Draw("COLZ");
	plot1->SetLineColor(kRed);
    	plot1->SetLineWidth(2);
	
	plot1->Draw("SAME");
	c->Update();
	
	system(Form("mkdir %d",n_channel));
	c->SaveAs(Form("%d/fprompt.pdf",n_channel));
	
	//--------------------------------
	//nessa segunda parte é para plotar o histograma de carga e fitar a carga
	hQ->Draw();
	double fit_min=1200e3;
	double fit_max=1400e3;
	
	double param[3];
	param[0]=3000;
	param[1]=(fit_min+fit_max)/2;
	param[2]=75e3;
	auto fit = new TF1("gaus_fit","gaus",fit_min,fit_max);
	fit->SetParameters(param);
	hQ->Fit("gaus_fit","0R");
	fit->SetRange(min_bin,max_bin);
	fit->Draw("SAME");
	
	fit->GetParameters(param);
	//legenda
	auto legend = new TLegend(0.6,0.7,1,0.9);
	legend->SetHeader("Legend","C"); // option "C" allows to center the header
	
	legend->AddEntry(hQ,"charge_histogram");
	legend->AddEntry(fit,Form("fit: -- AO: %0.2f -- mean: %0.2f -- std: %0.2f --",param[0],param[1],param[2]));
	legend->Draw("SAME");
	
	c->SaveAs(Form("%d/alpha.pdf",n_channel));
}
