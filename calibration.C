const int N_ENERGIES = 5;
const int N_SH_STARTS = 3;

std::array<int, N_ENERGIES> ENERGIES = {100,140,180,240,300};
std::array<int, N_ENERGIES> RUNS_CALIB = {100677,100633,100671,100648,100639};
std::array<int, N_ENERGIES> RUNS_TEST = {100631,100673,100672,100647,100645};

std::array<double, N_ENERGIES*3> X_MAX = {4500,5000,6500, 5500,7000,8500, 8000,8500,10000, 10000,11000,13000, 11000,12000,14000}; 
std::array<double, N_ENERGIES*3> X_MIN = {1000,1500,4000, 2000,3000,5000, 2000,3500,6000, 2500,4000,8500, 3000,5000,9500};

/*
std::array<double, N_ENERGIES*3> X_MAX = {6000,7000,8000, 8000,8500,10500, 11000,12000,13000, 13000,13500,15000, 13000,15000,16000}; 
std::array<double, N_ENERGIES*3> X_MIN = {500,500,1500, 1000,1500,3500, 500,1000,3000, 500,1000,4000, 1000,1500,5000};
*/
/*
std::array<double, N_ENERGIES*3> X_MAX = {4000,5000,6000, 6000,7500,9000, 9000,9500,11000, 10000,11000,13000, 11000,13000,14000}; 
std::array<double, N_ENERGIES*3> X_MIN = {1000,1500,3000, 2000,3000,5000, 1500,3500,6000, 1500,4000,7000, 3000,3000,9000};  

std::array<double, N_ENERGIES*3> X_MAX = {5000,6000,7000, 8000,8000,9000, 10000,11000,12000, 12000,12000,14000, 13500,15000,15000}; 
std::array<double, N_ENERGIES*3> X_MIN = {500,1000,2500, 1000,1000,3000, 1000,1000,4000, 1000,1000,4000, 1500,3000,6000};
*/
std::array<double, N_ENERGIES*3> A_MAX = {27,25,22, 21,21,20, 21,21,20, 21,21,25, 21,21,25};
std::array<double, N_ENERGIES*3> A_MIN = {3,3,5, 3,3,10, 3,3,10, 3,3,12, 3,3,15};
std::array<EColor, N_ENERGIES> COLORS = {kRed, kBlue, kGreen, kCyan, kYellow};

void styleGraph(TGraphErrors* &graph, const char* title, EColor color) {
  graph->SetTitle(title);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(21);
}

void styleMultigraph(TMultiGraph* &multigraph, const char* xlabel, const char* ylabel, double xmin, double xmax, double ymin, double ymax) {
  multigraph->GetXaxis()->SetTitle(xlabel);
  multigraph->GetYaxis()->SetTitle(ylabel);
  multigraph->GetXaxis()->SetLimits(xmin,xmax);
  multigraph->SetMinimum(ymin);
  multigraph->SetMaximum(ymax);
}

void openFiles(std::map<std::string, TFile*> &inFiles) {
  for (int i{0}; i<N_ENERGIES; ++i) {
    inFiles[Form("Calib_%i_GeV", ENERGIES[i])] = new TFile(Form("calib_alex/TB_outputRun_%i_ClusterSize_3.root",RUNS_CALIB[i]),"read");
    inFiles[Form("Test_%i_GeV", ENERGIES[i])] = new TFile(Form("calib_alex/TB_outputRun_%i_ClusterSize_3.root",RUNS_TEST[i]),"read");
  }
}

void definePlots(std::map<std::string, TH1*> &histos, std::map<std::string, TF1*> &functions, std::map<std::string, TGraphErrors*> &graphs, std::map<std::string, TMultiGraph*> &multigraphs) {

  // vs E
  multigraphs["SciFiQDC_vs_E"] = new TMultiGraph();
  styleMultigraph(multigraphs["SciFiQDC_vs_E"], "pion energy (GeV)", "mean SciFi QDC (a.u.)", 50, 350, 0, 5000);
  multigraphs["USQDC_vs_E"] = new TMultiGraph();
  styleMultigraph(multigraphs["USQDC_vs_E"], "pion energy (GeV)", "mean US QDC (a.u.)", 50, 350, 0, 14000);
  multigraphs["K_vs_E"] = new TMultiGraph();
  styleMultigraph(multigraphs["K_vs_E"], "pion energy (GeV)", "K (GeV/QDC)", 50, 350, 0.03, 0.15);
  multigraphs["Alpha_vs_E"] = new TMultiGraph();
  styleMultigraph(multigraphs["Alpha_vs_E"], "pion energy (GeV)", "Alpha (GeV/QDC)", 50, 350, -0.002, 0.025);
  multigraphs["OffsetE_vs_E"] = new TMultiGraph();
  styleMultigraph(multigraphs["OffsetE_vs_E"], "pion energy (GeV)", "E_reco - E (GeV)", 50, 350, -100, 100);
  multigraphs["OffsetErel_vs_E"] = new TMultiGraph();
  styleMultigraph(multigraphs["OffsetErel_vs_E"], "pion energy (GeV)", "(E_reco - E)/E %", 50, 350, -100, 100);
  multigraphs["Eres_vs_E"] = new TMultiGraph();
  styleMultigraph(multigraphs["Eres_vs_E"], "pion energy (GeV)", "Resolution %", 50, 350, 0, 50);
  
  // vs shStart
  multigraphs["SciFiQDC_vs_shStart"] = new TMultiGraph();
  styleMultigraph(multigraphs["SciFiQDC_vs_shStart"], "shower start", "mean SciFi QDC (a.u.)", 0.5, 4.5, 0, 5000);
  multigraphs["USQDC_vs_shStart"] = new TMultiGraph();
  styleMultigraph(multigraphs["USQDC_vs_shStart"], "shower start", "mean US QDC (a.u.)", 0.5, 4.5, 0, 14000);
  multigraphs["K_vs_shStart"] = new TMultiGraph();
  styleMultigraph(multigraphs["K_vs_shStart"], "shower start", "K (GeV/QDC)", 0.5, 4.5, 0.03, 0.15);
  multigraphs["Alpha_vs_shStart"] = new TMultiGraph();
  styleMultigraph(multigraphs["Alpha_vs_shStart"], "shower start", "Alpha (GeV/QDC)", 0.5, 4.5, -0.002, 0.025);
  multigraphs["OffsetE_vs_shStart"] = new TMultiGraph();
  styleMultigraph(multigraphs["OffsetE_vs_shStart"], "shower start", "E_reco - E (GeV)", 0.5, 4.5, -100, 100);
  multigraphs["OffsetErel_vs_shStart"] = new TMultiGraph();
  styleMultigraph(multigraphs["OffsetErel_vs_shStart"], "shower start", "(E_reco - E)/E %", 0.5, 4.5, -100, 100);
  multigraphs["Eres_vs_shStart"] = new TMultiGraph();
  styleMultigraph(multigraphs["Eres_vs_shStart"], "shower start", "Resolution %", 0.5, 4.5, 0, 50);

  for (int k{0}; k<N_SH_STARTS; ++k) {
    auto label_shStart = Form("start_in_SciFi%i",k+2);

    graphs[Form("SciFiQDC_vs_E_shStart%i",k+2)] = new TGraphErrors();
    styleGraph(graphs[Form("SciFiQDC_vs_E_shStart%i",k+2)], label_shStart, COLORS[k]);
    multigraphs["SciFiQDC_vs_E"]->Add(graphs[Form("SciFiQDC_vs_E_shStart%i",k+2)]);

    graphs[Form("USQDC_vs_E_shStart%i",k+2)] = new TGraphErrors();
    styleGraph(graphs[Form("USQDC_vs_E_shStart%i",k+2)], label_shStart, COLORS[k]);
    multigraphs["USQDC_vs_E"]->Add(graphs[Form("USQDC_vs_E_shStart%i",k+2)]);

    graphs[Form("K_vs_E_shStart%i",k+2)] = new TGraphErrors();
    styleGraph(graphs[Form("K_vs_E_shStart%i",k+2)], label_shStart, COLORS[k]);
    multigraphs["K_vs_E"]->Add(graphs[Form("K_vs_E_shStart%i",k+2)]);

    graphs[Form("Alpha_vs_E_shStart%i",k+2)] = new TGraphErrors();
    styleGraph(graphs[Form("Alpha_vs_E_shStart%i",k+2)], label_shStart, COLORS[k]);
    multigraphs["Alpha_vs_E"]->Add(graphs[Form("Alpha_vs_E_shStart%i",k+2)]);

    graphs[Form("OffsetE_vs_E_shStart%i",k+2)] = new TGraphErrors();
    styleGraph(graphs[Form("OffsetE_vs_E_shStart%i",k+2)], label_shStart, COLORS[k]);
    multigraphs["OffsetE_vs_E"]->Add(graphs[Form("OffsetE_vs_E_shStart%i",k+2)]);

    graphs[Form("OffsetErel_vs_E_shStart%i",k+2)] = new TGraphErrors();
    styleGraph(graphs[Form("OffsetErel_vs_E_shStart%i",k+2)], label_shStart, COLORS[k]);
    multigraphs["OffsetErel_vs_E"]->Add(graphs[Form("OffsetErel_vs_E_shStart%i",k+2)]);

    graphs[Form("Eres_vs_E_shStart%i",k+2)] = new TGraphErrors();
    styleGraph(graphs[Form("Eres_vs_E_shStart%i",k+2)], label_shStart, COLORS[k]);
    multigraphs["Eres_vs_E"]->Add(graphs[Form("Eres_vs_E_shStart%i",k+2)]);
  }

  for (int k{0}; k<N_ENERGIES; ++k) {
    auto label_En = Form("%i GeV",ENERGIES[k]);
    graphs[Form("SciFiQDC_vs_shStart_%i_GeV",ENERGIES[k])] = new TGraphErrors();
    styleGraph(graphs[Form("SciFiQDC_vs_shStart_%i_GeV",ENERGIES[k])], label_En, COLORS[k]);
    multigraphs["SciFiQDC_vs_shStart"]->Add(graphs[Form("SciFiQDC_vs_shStart_%i_GeV",ENERGIES[k])]);

    graphs[Form("USQDC_vs_shStart_%i_GeV",ENERGIES[k])] = new TGraphErrors();
    styleGraph(graphs[Form("USQDC_vs_shStart_%i_GeV",ENERGIES[k])], label_En, COLORS[k]);
    multigraphs["USQDC_vs_shStart"]->Add(graphs[Form("USQDC_vs_shStart_%i_GeV",ENERGIES[k])]);

    graphs[Form("K_vs_shStart_%i_GeV",ENERGIES[k])] = new TGraphErrors();
    styleGraph(graphs[Form("K_vs_shStart_%i_GeV",ENERGIES[k])], label_En, COLORS[k]);
    multigraphs["K_vs_shStart"]->Add(graphs[Form("K_vs_shStart_%i_GeV",ENERGIES[k])]);

    graphs[Form("Alpha_vs_shStart_%i_GeV",ENERGIES[k])] = new TGraphErrors();
    styleGraph(graphs[Form("Alpha_vs_shStart_%i_GeV",ENERGIES[k])], label_En, COLORS[k]);
    multigraphs["Alpha_vs_shStart"]->Add(graphs[Form("Alpha_vs_shStart_%i_GeV",ENERGIES[k])]);

    graphs[Form("OffsetE_vs_shStart_%i_GeV",ENERGIES[k])] = new TGraphErrors();
    styleGraph(graphs[Form("OffsetE_vs_shStart_%i_GeV",ENERGIES[k])], label_En, COLORS[k]);
    multigraphs["OffsetE_vs_shStart"]->Add(graphs[Form("OffsetE_vs_shStart_%i_GeV",ENERGIES[k])]);

    graphs[Form("OffsetErel_vs_shStart_%i_GeV",ENERGIES[k])] = new TGraphErrors();
    styleGraph(graphs[Form("OffsetErel_vs_shStart_%i_GeV",ENERGIES[k])], label_En, COLORS[k]);
    multigraphs["OffsetErel_vs_shStart"]->Add(graphs[Form("OffsetErel_vs_shStart_%i_GeV",ENERGIES[k])]);

    graphs[Form("Eres_vs_shStart_%i_GeV",ENERGIES[k])] = new TGraphErrors();
    styleGraph(graphs[Form("Eres_vs_shStart_%i_GeV",ENERGIES[k])], label_En, COLORS[k]);
    multigraphs["Eres_vs_shStart"]->Add(graphs[Form("Eres_vs_shStart_%i_GeV",ENERGIES[k])]);

    histos[Form("Scifi/E_reco vs starting wall %i GeV",ENERGIES[k])] = new TH2D(Form("Scifi/E_reco vs starting wall %i GeV",ENERGIES[k]),Form("Scifi/E_reco vs starting wall %i GeV; starting wall; SciFi/E_reco (a.u.)", ENERGIES[k]), 3, 0.5, 3.5, 200, -1, 1);
  }
    
  for (int m{0}; m<N_ENERGIES; ++m) {
    for (int l{0}; l<N_SH_STARTS; ++l) {
      histos[Form("Alpha_%i_GeV_shStart%i",ENERGIES[m],l+2)] = new TH1D(Form("Alpha_%i_GeV_shStart%i",ENERGIES[m],l+2),Form("Alpha_%i_GeV_shStart%i; Alpha (GeV/QDC); entries",ENERGIES[m],l+2), 70, -0.02, 0.04);
      histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[m],l+2)] = new TH1D(Form("Ereco_%i_GeV_shStart%i",ENERGIES[m],l+2),Form("Ereco_%i_GeV_shStart%i; Reconstructed energy (GeV); entries",ENERGIES[m],l+2), 100, 0, 500);
      histos[Form("Ereco_%i_GeV_shStart%i_low",ENERGIES[m],l+2)] = new TH1D(Form("Ereco_%i_GeV_shStart%i_low",ENERGIES[m],l+2),Form("Ereco_%i_GeV_shStart%i_low; Reconstructed energy (GeV); entries",ENERGIES[m],l+2), 100, 0, 500);
      histos[Form("Ereco_%i_GeV_shStart%i_high",ENERGIES[m],l+2)] = new TH1D(Form("Ereco_%i_GeV_shStart%i_high",ENERGIES[m],l+2),Form("Ereco_%i_GeV_shStart%i_high; Reconstructed energy (GeV); entries",ENERGIES[m],l+2), 100, 0, 500);
      histos[Form("high-low_%i_GeV_shStart%i",ENERGIES[m],l+2)] = new TH1D(Form("high-low_%i_GeV_shStart%i",ENERGIES[m],l+2),Form("high-low_%i_GeV_shStart%i; Energy gap (GeV); entries",ENERGIES[m],l+2), 1000, -100, 100);
      functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[m],l+2)] = new TF1(Form("f_%i_GeV_shStart%i", ENERGIES[m], l+2), "[0]*x + [1]");
      functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[m],l+2)]->SetParameters(-1/4, 4000);
      //functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[m],l+2)]->SetRange(X_MIN[3*m+l], X_MAX[3*m+l]);
      functions[Form("fitGaus_Alpha_%i_GeV_shStart%i",ENERGIES[m],l+2)] = new TF1("gaus", "gaus");
      functions[Form("fitGaus_Alpha_%i_GeV_shStart%i",ENERGIES[m],l+2)]->SetParameters(1800, 12*0.001, 4*0.001);
      functions[Form("fitGaus_Alpha_%i_GeV_shStart%i",ENERGIES[m],l+2)]->SetRange(A_MIN[3*m+l]*0.001, A_MAX[3*m+l]*0.001);
    }
  }

}

void fillQDCPLots(std::map<std::string, TFile*> &inFiles, std::map<std::string, TH1*> &histos, std::map<std::string, TGraphErrors*> &graphs, std::map<std::string, TProfile*> &profiles) {
  for (int i{0}; i<N_ENERGIES; ++i) {
    histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV",ENERGIES[i])] = (TH2D*)((TH2D*)inFiles[Form("Calib_%i_GeV", ENERGIES[i])]->Get("Cut_QDCUS_vs_QDCScifi"))->Clone(Form("QDCUS_vs_QDCScifi_Calib_%i_GeV", ENERGIES[i]));
    histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV",ENERGIES[i])] = (TH2D*)((TH2D*)inFiles[Form("Test_%i_GeV", ENERGIES[i])]->Get("Cut_QDCUS_vs_QDCScifi"))->Clone(Form("QDCUS_vs_QDCScifi_Test_%i_GeV", ENERGIES[i]));
    for (int j{0}; j<N_SH_STARTS; ++j) {
      histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH2D*)((TH2D*)inFiles[Form("Calib_%i_GeV", ENERGIES[i])]->Get(Form("Cut_QDCUS_vs_QDCScifi_ShStart_st%i",j+2)))->Clone(Form("QDCUS_vs_QDCScifi_Calib_%i_GeV_st%i", ENERGIES[i], j+2));
      histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH2D*)((TH2D*)inFiles[Form("Test_%i_GeV", ENERGIES[i])]->Get(Form("Cut_QDCUS_vs_QDCScifi_ShStart_st%i",j+2)))->Clone(Form("QDCUS_vs_QDCScifi_Test_%i_GeV_st%i", ENERGIES[i], j+2));

      profiles[Form("Calib_profile_%i_GeV_shStart%i",ENERGIES[i],j+2)] =  ((TH2D*)histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)])->ProfileX(Form("Calib_profile_%i_GeV_shStart%i",ENERGIES[i],j+2));

      double scifiQDC = ((TH1F*)inFiles[Form("Calib_%i_GeV", ENERGIES[i])]->Get(Form("Cut_Shower_SciFi_QDC_shStart%i",j+2)))->GetMean();
      double scifiQDCerr = ((TH1F*)inFiles[Form("Calib_%i_GeV", ENERGIES[i])]->Get(Form("Cut_Shower_SciFi_QDC_shStart%i",j+2)))->GetStdDev();

      TH2D* h = (TH2D*)inFiles[Form("Calib_%i_GeV", ENERGIES[i])]->Get(Form("Cut_QDCUS_vs_QDCScifi_ShStart_st%i",j+2));
      int ny = h->GetYaxis()->GetNbins();

      double usQDC = ((TH1D*)h->ProjectionX("h1D", 0, ny))->GetMean();
      double usQDCerr = ((TH1D*)h->ProjectionX("h1D", 0, ny))->GetStdDev();

      graphs[Form("SciFiQDC_vs_E_shStart%i",j+2)]->SetPoint(i, ENERGIES[i] + 5*j, scifiQDC);
      graphs[Form("SciFiQDC_vs_E_shStart%i",j+2)]->SetPointError(i, 0, scifiQDCerr);
      graphs[Form("SciFiQDC_vs_shStart_%i_GeV",ENERGIES[i])]->SetPoint(j, j+2+0.1*i, scifiQDC);
      graphs[Form("SciFiQDC_vs_shStart_%i_GeV",ENERGIES[i])]->SetPointError(j, 0, scifiQDCerr);

      graphs[Form("USQDC_vs_E_shStart%i",j+2)]->SetPoint(i, ENERGIES[i] + 5*j, usQDC);
      graphs[Form("USQDC_vs_E_shStart%i",j+2)]->SetPointError(i, 0, usQDCerr);
      graphs[Form("USQDC_vs_shStart_%i_GeV",ENERGIES[i])]->SetPoint(j, j+2+0.1*i, usQDC);
      graphs[Form("USQDC_vs_shStart_%i_GeV",ENERGIES[i])]->SetPointError(j, 0, usQDCerr);
    }
  }
}

double computeKfromFit(std::map<std::string, TFile*> &inFiles, std::map<std::string, TH1*> &histos, std::map<std::string, TF1*> &functions, std::map<std::string, TGraphErrors*> &graphs) {
  double sum_k{0};
  for (int i{0}; i<N_ENERGIES; ++i) {
    for (int j{0}; j<N_SH_STARTS; ++j) {     
      // fit previously filled scatter plots

      functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->SetRange(X_MIN[3*i+j], X_MAX[3*i+j]);
      histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->Fit(functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)],"RQ");
      double k = ENERGIES[i]/functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParameter(1);
      double kErr = k*functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParError(1)/functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParameter(1);
      std::cout<<"k:\t"<<k<<std::endl;
      std::cout<<"errk:\t"<<kErr<<std::endl;
      graphs[Form("K_vs_E_shStart%i",j+2)]->SetPoint(i,ENERGIES[i] + 5*j,k);
      graphs[Form("K_vs_E_shStart%i",j+2)]->SetPointError(i,0,kErr);
      graphs[Form("K_vs_shStart_%i_GeV",ENERGIES[i])]->SetPoint(j,j+2+0.1*i,k);
      graphs[Form("K_vs_shStart_%i_GeV",ENERGIES[i])]->SetPointError(j,0,kErr);
      if (j<N_SH_STARTS-1) {
        sum_k += k;
      }
    }
  }
  return sum_k/((N_SH_STARTS-1)*N_ENERGIES);      
}

double computeAlphafromFit(std::map<std::string, TFile*> &inFiles, std::map<std::string, TH1*> &histos, std::map<std::string, TF1*> &functions, std::map<std::string, TGraphErrors*> &graphs) {
  double sum_alpha{0};
  for (int i{0}; i<N_ENERGIES; ++i) {
    for (int j{0}; j<N_SH_STARTS; ++j) {    
      // fit previously filled scatter plots
      functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->SetRange(X_MIN[3*i+j], X_MAX[3*i+j]);
      histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->Fit(functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)],"RQ");
      double alpha = -ENERGIES[i]*functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParameter(0)/functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParameter(1);
      
      double alphaErr = alpha*(functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParError(0)/functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParameter(0) + functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParError(1)/functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParameter(1));
      std::cout<<"alpha:\t"<<alpha<<std::endl;
      std::cout<<"erralpha:\t"<<alphaErr<<std::endl;
      graphs[Form("Alpha_vs_E_shStart%i",j+2)]->SetPoint(i,ENERGIES[i] + 5*j,alpha);
      graphs[Form("Alpha_vs_E_shStart%i",j+2)]->SetPointError(i,0,alphaErr);
      graphs[Form("Alpha_vs_shStart_%i_GeV",ENERGIES[i])]->SetPoint(j,j+2+0.1*i,alpha);
      graphs[Form("Alpha_vs_shStart_%i_GeV",ENERGIES[i])]->SetPointError(j,0,alphaErr);
      if (j<N_SH_STARTS-1) {
        sum_alpha += alpha;
      }
    }
  }
  return sum_alpha/((N_SH_STARTS-1)*N_ENERGIES);      
}

double computeAlpha(std::map<std::string, TFile*> &inFiles, std::map<std::string, TH1*> &histos, std::map<std::string, TF1*> &functions, std::map<std::string, TGraphErrors*> &graphs, const double k) {
  double sum_alpha{0};
  for (int i{0}; i<N_ENERGIES; ++i) {
    for (int j{0}; j<N_SH_STARTS; ++j) {  
      int upper_x_bin = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->FindBin(X_MAX[3*i+j]);
      int lower_x_bin = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->FindBin(X_MIN[3*i+j]);
      int upper_y_bin = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetYaxis()->GetNbins();
      for (int binX = lower_x_bin; binX <= upper_x_bin; ++binX) {
        for (int binY = 1; binY <= upper_y_bin; ++binY) {
            // Access the bin content
            double binContent = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetBinContent(binX, binY);
            if (binContent>0) {
              double scifiQDC = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetYaxis()->GetBinCenter(binY);
              double usQDC = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->GetBinCenter(binX);
              double a = (ENERGIES[i] - k * scifiQDC)/usQDC;
              histos[Form("Alpha_%i_GeV_shStart%i",ENERGIES[i],j+2)]->Fill(a,binContent);             
            }
        }
      }   
      // fit obtained alpha distibution
      histos[Form("Alpha_%i_GeV_shStart%i",ENERGIES[i],j+2)]->Fit(functions[Form("fitGaus_Alpha_%i_GeV_shStart%i",ENERGIES[i],j+2)],"RQ");
      double alpha = functions[Form("fitGaus_Alpha_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParameter(1);
      double alphaErr = functions[Form("fitGaus_Alpha_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParameter(2); // Suggestion: take the error on the mean (1st parameter) as the error instead
      std::cout<<"alpha:\t"<<alpha<<std::endl;
      std::cout<<"alphaErr:\t"<<alphaErr<<std::endl;
      graphs[Form("Alpha_vs_E_shStart%i",j+2)]->SetPoint(i,ENERGIES[i] + 5*j,alpha);
      graphs[Form("Alpha_vs_E_shStart%i",j+2)]->SetPointError(i,0,alphaErr);
      graphs[Form("Alpha_vs_shStart_%i_GeV",ENERGIES[i])]->SetPoint(j,j+2+0.1*i,alpha);
      graphs[Form("Alpha_vs_shStart_%i_GeV",ENERGIES[i])]->SetPointError(j,0,alphaErr);
      if (j<N_SH_STARTS-1) {
        sum_alpha += alpha;
      }
    }
  }
  return sum_alpha/((N_SH_STARTS-1)*N_ENERGIES);      
}

void pca(std::map<std::string, TFile*> &inFiles, std::map<std::string, TH1*> &histos, std::map<std::string, TGraphErrors*> &graphs, double &k_mean, double &alpha_mean, std::map<std::string, TLine*> &pca_lines, std::vector<double> &ks, std::vector<double> &alphas){
  double* data = new double[2];
  double alpha_sum = 0;
  double k_sum = 0;
  for (int i{0}; i<N_ENERGIES; ++i) {
    for (int j{0}; j<N_SH_STARTS; ++j) { 
      TPrincipal* principal = new TPrincipal(2,"D"); 
      int upper_x_bin = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->FindBin(X_MAX[3*i+j]);
      int lower_x_bin = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->FindBin(X_MIN[3*i+j]);
      int upper_y_bin = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetYaxis()->GetNbins();
      for (int binX = lower_x_bin; binX <= upper_x_bin; ++binX) {
        for (int binY = 1; binY <= upper_y_bin; ++binY) {
          // Access the bin content
          double binContent = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetBinContent(binX, binY);
          if (binContent>0) {
            double scifiQDC = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetYaxis()->GetBinCenter(binY);
            double usQDC = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->GetBinCenter(binX);
            data[0] = usQDC;
            data[1] = scifiQDC;
            //std::cout<<data[0]<<"\t"<<data[1]<<"\n";
            for (int m{0}; m<binContent; ++m) {
              principal->AddRow(data);
            }            
          }
        }
      }   
      // make pca
      //std::cout<<"\n\n\nEnergy:\t"<<ENERGIES[i]<<"\t station:\t"<< j+1 << "\n";
      principal->MakePrincipals();
      // principal->Print();
      auto eigenV = principal->GetEigenVectors();
      // eigenV->Print();
      auto meanV = principal->GetMeanValues();
      // meanV->Print();
      auto sigmaV = principal->GetSigmas();
      // principal->MakeHistograms();
      // auto hist_list = principal->GetHistograms();
      // histos[Form("pca_s_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH1D*)((TH1D*)hist_list->FindObject("pca_s"))->Clone(Form("pca_s_%i_GeV_shStart%i",ENERGIES[i],j+2));
      // histos[Form("pca_e_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH1D*)((TH1D*)hist_list->FindObject("pca_e"))->Clone(Form("pca_e_%i_GeV_shStart%i",ENERGIES[i],j+2));
      // histos[Form("pca_x000_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH1D*)((TH1D*)hist_list->FindObject("pca_x000"))->Clone(Form("pca_x000_%i_GeV_shStart%i",ENERGIES[i],j+2));
      // histos[Form("pca_d000_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH1D*)((TH1D*)hist_list->FindObject("pca_d000"))->Clone(Form("pca_d000_%i_GeV_shStart%i",ENERGIES[i],j+2));
      // histos[Form("pca_p000_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH1D*)((TH1D*)hist_list->FindObject("pca_p000"))->Clone(Form("pca_p000_%i_GeV_shStart%i",ENERGIES[i],j+2));
      // histos[Form("pca_x001_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH1D*)((TH1D*)hist_list->FindObject("pca_x001"))->Clone(Form("pca_x001_%i_GeV_shStart%i",ENERGIES[i],j+2));
      // histos[Form("pca_d001_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH1D*)((TH1D*)hist_list->FindObject("pca_d001"))->Clone(Form("pca_d001_%i_GeV_shStart%i",ENERGIES[i],j+2));
      // histos[Form("pca_p001_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH1D*)((TH1D*)hist_list->FindObject("pca_p001"))->Clone(Form("pca_p001_%i_GeV_shStart%i",ENERGIES[i],j+2));
      double m = (*eigenV)[0][1]/(*eigenV)[0][0];
      double q = (*meanV)[1]-m*(*meanV)[0];
      pca_lines[Form("Line_%i_GeV_shStart%i",ENERGIES[i],j+2)] = new TLine(X_MIN[3*i+j], m*X_MIN[3*i+j]+q, X_MAX[3*i+j], m*X_MAX[3*i+j]+q);
      double k = ENERGIES[i]/q;
      double alpha = -m*k;

      double dk_dy = -ENERGIES[i] / pow(((*meanV)[0] - m * (*meanV)[1]), 2);
      double dk_dx = (-ENERGIES[i] * m) / pow(((*meanV)[0] - m * (*meanV)[1]), 2);

      double kErr = std::sqrt(std::pow(dk_dy * (*sigmaV)[1], 2)+ std::pow(dk_dx * (*sigmaV)[0], 2));
      double alphaErr = std::abs(m*kErr);
      //principal->GetSigmas()->Print();
      //principal->MakeHistograms();
      delete principal;
      
      // double alphaErr = (*std::max_element(alpha_temp.begin(),alpha_temp.end()) - *std::min_element(alpha_temp.begin(),alpha_temp.end()))/2;
      // double kErr = (*std::max_element(k_temp.begin(),k_temp.end()) - *std::min_element(k_temp.begin(),k_temp.end()))/2;
      graphs[Form("Alpha_vs_E_shStart%i",j+2)]->SetPoint(i,ENERGIES[i] + 5*j,alpha);
      graphs[Form("Alpha_vs_E_shStart%i",j+2)]->SetPointError(i,0,alphaErr);
      graphs[Form("Alpha_vs_shStart_%i_GeV",ENERGIES[i])]->SetPoint(j,j+2+0.1*i,alpha);
      graphs[Form("Alpha_vs_shStart_%i_GeV",ENERGIES[i])]->SetPointError(j,0,alphaErr);

      graphs[Form("K_vs_E_shStart%i",j+2)]->SetPoint(i,ENERGIES[i] + 5*j,k);
      graphs[Form("K_vs_E_shStart%i",j+2)]->SetPointError(i,0,kErr);
      graphs[Form("K_vs_shStart_%i_GeV",ENERGIES[i])]->SetPoint(j,j+2+0.1*i,k);
      graphs[Form("K_vs_shStart_%i_GeV",ENERGIES[i])]->SetPointError(j,0,kErr);
      if (j<N_SH_STARTS-1) {
        ks.push_back(k);
        alphas.push_back(alpha);
        alpha_sum += alpha;
        k_sum += k;
      } 
    }
  }
  delete [] data;
  alpha_mean = alpha_sum/((N_SH_STARTS-1)*N_ENERGIES);
  k_mean = k_sum/((N_SH_STARTS-1)*N_ENERGIES);
}

void testEnergyResolution(std::map<std::string, TFile*> &inFiles, std::map<std::string, TH1*> &histos, std::map<std::string, TGraphErrors*> &graphs, const double k, const double alpha, const std::vector<double> ks, const std::vector<double> alphas, std::map<std::string, TF1*> &functions) {
  // const double k_min = *std::min_element(ks.begin(), ks.end());
  // const double k_max = *std::max_element(ks.begin(), ks.end());
  // const double alpha_min = *std::min_element(alphas.begin(), alphas.end());
  // const double alpha_max = *std::max_element(alphas.begin(), alphas.end());

  const double k_min = ks[0];
  const double k_max = ks[8];
  const double alpha_min = alphas[0];
  const double alpha_max = alphas[8];

  for (int i{0}; i<N_ENERGIES; ++i) {
    for (int j{0}; j<N_SH_STARTS; ++j) { 
      // int upper_x_bin = histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->GetNbins();
      // int upper_y_bin = histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetYaxis()->GetNbins();
      int upper_x_bin = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->FindBin(X_MAX[3*i+j]);
      int lower_x_bin = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->FindBin(X_MIN[3*i+j]);
      int upper_y_bin = histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetYaxis()->GetNbins();
      for (int binX = lower_x_bin; binX <= upper_x_bin; ++binX) {
        for (int binY = 1; binY <= upper_y_bin; ++binY) {
          // Access the bin content
          double binContent = histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetBinContent(binX, binY);
          if (binContent>0) {
            double scifiQDC = histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetYaxis()->GetBinCenter(binY);
            double usQDC = histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->GetBinCenter(binX);
            double recEn = k*scifiQDC + alpha*usQDC;
            double recEn_min = k_min*scifiQDC + alpha_min*usQDC;
            double recEn_max = k_max*scifiQDC + alpha_max*usQDC;
            histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[i],j+2)]->Fill(recEn,binContent);
            histos[Form("Ereco_%i_GeV_shStart%i_low",ENERGIES[i],j+2)]->Fill(recEn_min,binContent);
            histos[Form("Ereco_%i_GeV_shStart%i_high",ENERGIES[i],j+2)]->Fill(recEn_max,binContent);
            histos[Form("high-low_%i_GeV_shStart%i",ENERGIES[i],j+2)]->Fill(recEn_max - recEn_min,binContent);
            if (recEn>0) {
              histos[Form("Scifi/E_reco vs starting wall %i GeV",ENERGIES[i])]->Fill(j+1, scifiQDC*k/recEn);
            }            
          }
        }
      }
      double h_max = histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->GetBinCenter(histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetMaximumBin());
      double h_stddev = histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetStdDev();

      functions[Form("Gaus_%i_GeV_shStart%i",ENERGIES[i],j+2)] = new TF1(Form("Gaus_%i_GeV_shStart%i",ENERGIES[i],j+2), "gaus", h_max-1.5*h_stddev, h_max+1.5*h_stddev);
      histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[i],j+2)]->Fit(Form("Gaus_%i_GeV_shStart%i",ENERGIES[i],j+2), "RQ");
      double Ereco = functions[Form("Gaus_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParameter(1);
      double ErecoErr = functions[Form("Gaus_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetParameter(2);

      graphs[Form("OffsetE_vs_E_shStart%i",j+2)]->SetPoint(i,ENERGIES[i] + 5*j,Ereco - ENERGIES[i]);
      graphs[Form("OffsetE_vs_E_shStart%i",j+2)]->SetPointError(i,0,ErecoErr);
      graphs[Form("OffsetE_vs_shStart_%i_GeV",ENERGIES[i])]->SetPoint(j,j+2+0.1*i,Ereco - ENERGIES[i]);
      graphs[Form("OffsetE_vs_shStart_%i_GeV",ENERGIES[i])]->SetPointError(j,0,ErecoErr);

      graphs[Form("OffsetErel_vs_E_shStart%i",j+2)]->SetPoint(i,ENERGIES[i] + 5*j,(Ereco - ENERGIES[i])/ENERGIES[i]*100);
      graphs[Form("OffsetErel_vs_E_shStart%i",j+2)]->SetPointError(i,0,ErecoErr/ENERGIES[i]*100);
      graphs[Form("OffsetErel_vs_shStart_%i_GeV",ENERGIES[i])]->SetPoint(j,j+2+0.1*i,(Ereco - ENERGIES[i])/ENERGIES[i]*100);
      graphs[Form("OffsetErel_vs_shStart_%i_GeV",ENERGIES[i])]->SetPointError(j,0,ErecoErr/ENERGIES[i]*100);

      graphs[Form("Eres_vs_E_shStart%i",j+2)]->SetPoint(i,ENERGIES[i] + 5*j,ErecoErr/ENERGIES[i]*100);
      graphs[Form("Eres_vs_shStart_%i_GeV",ENERGIES[i])]->SetPoint(j,j+2+0.1*i,ErecoErr/ENERGIES[i]*100);
    }
  }
}

void drawOnCanvases(std::map<std::string, TH1*> &histos, std::map<std::string, TF1*> &functions, std::map<std::string, TMultiGraph*> &multigraphs, std::map<std::string, TCanvas*> &canvases, std::map<std::string, TLine*> &pca_lines) {
  for (auto const& mg : multigraphs) {
    canvases[mg.first] = new TCanvas(mg.first.c_str(),mg.first.c_str(),1980,1020);
    canvases[mg.first]->cd();
    mg.second->Draw("ape");
    canvases[mg.first]->BuildLegend();
  }

  for (int k{0}; k<N_SH_STARTS; ++k) {
    for (int m{0}; m<N_ENERGIES; ++m) {
      canvases[Form("GausFit_%i_GeV_shStart%i",ENERGIES[m],k+2)] = new TCanvas(Form("GausFit_%i_GeV_shStart%i",ENERGIES[m],k+2),Form("GausFit_%i_GeV_shStart%i",ENERGIES[m],k+2),1980,1020);
      canvases[Form("GausFit_%i_GeV_shStart%i",ENERGIES[m],k+2)]->cd();
      gStyle->SetOptFit(1);
      histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[m],k+2)]->Draw("hist");
      functions[Form("Gaus_%i_GeV_shStart%i",ENERGIES[m],k+2)]->Draw("same");
      // TPaveStats *fitStats = (TPaveStats*)histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[m],k+2)]->FindObject("stats");
      // fitStats->Draw("same");
    }
  }

  for (int k{0}; k<N_SH_STARTS; ++k) {
    canvases[Form("Ereco_shStart%i",k+2)] = new TCanvas(Form("Ereco_shStart%i",k+2),Form("Ereco_shStart%i",k+2),1980,1020);
    canvases[Form("Ereco_shStart%i",k+2)]->cd();
    for (int m{0}; m<N_ENERGIES; ++m) {
      histos[Form("Clone_Ereco_%i_GeV_shStart%i",ENERGIES[m],k+2)] = (TH1D*)((TH1D*)histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[m],k+2)]->Clone(Form("Clone_Ereco_%i_GeV_shStart%i",ENERGIES[m],k+2)));
      histos[Form("Clone_Ereco_%i_GeV_shStart%i",ENERGIES[m],k+2)]->SetLineColor(COLORS[m]);
      histos[Form("Clone_Ereco_%i_GeV_shStart%i",ENERGIES[m],k+2)]->SetStats(0);
      histos[Form("Clone_Ereco_%i_GeV_shStart%i",ENERGIES[m],k+2)]->Scale(1.0 / histos[Form("Clone_Ereco_%i_GeV_shStart%i",ENERGIES[m],k+2)]->Integral());
      if (m==0) {
        histos[Form("Clone_Ereco_%i_GeV_shStart%i",ENERGIES[m],k+2)]->Draw("hist");
      }
      else {
        histos[Form("Clone_Ereco_%i_GeV_shStart%i",ENERGIES[m],k+2)]->Draw("hist same");
      }
    }
    canvases[Form("Ereco_shStart%i",k+2)]->BuildLegend();
  }

  for (int k{0}; k<N_SH_STARTS; ++k) {
    for (int m{0}; m<N_ENERGIES; ++m) {
      canvases[Form("PCA_%i_GeV_shStart%i",ENERGIES[m],k+2)] = new TCanvas(Form("PCA_%i_GeV_shStart%i",ENERGIES[m],k+2),Form("PCA_%i_GeV_shStart%i",ENERGIES[m],k+2),1980,1020);
      canvases[Form("PCA_%i_GeV_shStart%i",ENERGIES[m],k+2)]->cd();
      histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[m],k+2)]->Draw("COLZ");
      pca_lines[Form("Line_%i_GeV_shStart%i",ENERGIES[m],k+2)]->Draw();
    }
  }
}

void calibration() {
  std::map<std::string, TFile*> inFiles;
  std::map<std::string, TH1*> histos;
  std::map<std::string, TF1*> functions;
  std::map<std::string, TGraphErrors*> graphs;
  std::map<std::string, TMultiGraph*> multigraphs;
  std::map<std::string, TProfile*> profiles;
  std::map<std::string, TLine*> pca_lines;
  std::map<std::string, TCanvas*> canvases;
  double k{0};
  double alpha{0};
  std::vector<double> ks;
  std::vector<double> alphas;


  //TFile* outputFile = new TFile("calib_out_alex1/calibrations_pca.root","recreate");
  TFile* outputFile = new TFile("calibrations_pca.root","recreate");
  definePlots(histos, functions, graphs, multigraphs);
  openFiles(inFiles);

  fillQDCPLots(inFiles, histos, graphs, profiles);
  pca(inFiles, histos, graphs, k, alpha, pca_lines, ks, alphas);
  //k = computeKfromFit(inFiles, histos, functions, graphs);
  //alpha = computeAlphafromFit(inFiles, histos, functions, graphs);
  std::cout<<"mean K:\t"<<k<<std::endl;
  std::cout<<"mean Alpha:\t"<<alpha<<std::endl;
  std::cout<<"min K:\t"<<*std::min_element(ks.begin(), ks.end())<<std::endl;
  std::cout<<"max K:\t"<<*std::max_element(ks.begin(), ks.end())<<std::endl;
  std::cout<<"min Alpha:\t"<<*std::min_element(alphas.begin(), alphas.end())<<std::endl;
  std::cout<<"max Alpha:\t"<<*std::max_element(alphas.begin(), alphas.end())<<std::endl;
  testEnergyResolution(inFiles, histos, graphs, k, alpha, ks, alphas, functions);
  drawOnCanvases(histos, functions, multigraphs, canvases, pca_lines);

  outputFile->cd();

  for (auto const& h : histos) {
    h.second->Write();
  }
  for (auto const& c : canvases) {
    c.second->Write();
  }

  for (auto const& p : profiles) {
    p.second->Write();
  }
  
  //outputFile->Write();
  outputFile->Close();

}