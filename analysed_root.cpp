typedef struct
{
	double charge_1;
	double charge_2;
	double charge_t;
	double charge_p;
	double amplitude_1;
	double amplitude_2;
	double amplitude_t;
	double amplitude_p;
	double baseline_1;
	double baseline_2;
	double noise_1;
	double noise_2;
	double fprompt_1;
	double fprompt_2;
	double fprompt_t;
	long long n_event;
}event;


void analysed_root(string filename)
{
	TFile *f1 = new TFile(filename.c_str(),"READ"); //abre o .root, apenas le
	TTree *t1 = (TTree*)f1->Get("TChan"); //pega a ttree


	double baseline;
	long long index;
	double noise;	
	double roih;
	double roiq;
	long long nevent;
	double fprompt;
	
	t1->SetBranchAddress("baseline",&baseline);//aponta a variavel criada, para a localização correta na ttree
	t1->SetBranchAddress("index",&index);//aponta a variavel criada, para a localização correta na ttree
	t1->SetBranchAddress("noise",&noise);//aponta a variavel criada, para a localização correta na ttree
	t1->SetBranchAddress("roih",&roih);//aponta a variavel criada, para a localização correta na ttree
	t1->SetBranchAddress("roiq",&roiq);//aponta a variavel criada, para a localização correta na ttree
	t1->SetBranchAddress("nevent",&nevent);//aponta a variavel criada, para a localização correta na ttree
	t1->SetBranchAddress("fprompt",&fprompt);//aponta a variavel criada, para a localização correta na ttree
	
	
	Int_t entries=(Int_t)t1->GetEntries();
	
	//cria novo .root
	string new_name = filename.substr(0,filename.find(".root"))+"_analysed.root";
	TFile *f2 = new TFile(new_name.c_str(),"RECREATE");

	TTree* t2 = new TTree("t2","Data");
	event my_data;
	TBranch* tb= t2->Branch("data",&my_data,"charge_1/D:charge_2/D:charge_t/D:charge_p/D:amplitude_1/D:amplitude_2/D:amplitude_t/D:amplitude_p/D:baseline_1/D:baseline_2/D:noise_1/D:noise_2/D:fprompt_1/D:fprompt_2/D:fprompt_t/D:nevent/L");
	
	
	
	t1->GetEntry(0); //pega a primeira entrada
	long long old_nevent=nevent;
	
	
	//varre as entradas
	for(int i=0;i<entries;i++)
	{
		t1->GetEntry(i); //pega a i_esima entrada
		if(index==1)
		{
			my_data.charge_1=roiq;
			my_data.amplitude_1=roih;
			my_data.fprompt_1=fprompt;
			my_data.noise_1=noise;
			my_data.baseline_1=baseline;
		}
		else if (index==3)
		{
			my_data.charge_2=roiq;
			my_data.amplitude_2=roih;
			my_data.fprompt_2=fprompt;
			my_data.noise_2=noise;
			my_data.baseline_2=baseline;
			
			my_data.charge_t=my_data.charge_1+my_data.charge_2;
			my_data.fprompt_t=my_data.fprompt_1+my_data.fprompt_2;
			if(my_data.amplitude_1>my_data.amplitude_2)
				my_data.amplitude_t=my_data.amplitude_1;
			else
				my_data.amplitude_t=my_data.amplitude_2;
			
		}
		else if (index == 5)
		{
			my_data.charge_p=roiq;
			my_data.amplitude_p=roih;
			my_data.n_event=nevent;
			t2->Fill();
		}
		
		
	
	}
	
	f1->Close();
	f2->WriteObject(t2,"t2","TObject::kOverwrite");
	f2->Close();

}
