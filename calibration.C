const int N_ENERGIES = 5;
const int N_SH_STARTS = 3;
std::array<int, N_ENERGIES> ENERGIES = {100,140,180,240,300};
std::array<int, N_ENERGIES> RUNS_CALIB = {100677,100633,100671,100648,100639};
std::array<int, N_ENERGIES> RUNS_TEST = {100631,100673,100672,100647,100645};
std::array<double, N_ENERGIES*3> X_MAX = {5000,6000,7000, 8000,8000,9000, 10000,11000,12000, 12000,12000,14000, 13500,15000,15000}; 
std::array<double, N_ENERGIES*3> X_MIN = {500,1000,2500, 1000,1000,3000, 1000,1000,4000, 1000,1000,4000, 1500,3000,6000};  
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
    inFiles[Form("Calib_%i_GeV", ENERGIES[i])] = new TFile(Form("calib_collab/TB_outputRun_%i.root",RUNS_CALIB[i]),"read");
    inFiles[Form("Test_%i_GeV", ENERGIES[i])] = new TFile(Form("calib_collab/TB_outputRun_%i.root",RUNS_TEST[i]),"read");
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
      functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[m],l+2)] = new TF1(Form("f_%i_GeV_shStart%i", ENERGIES[m], l+2), "[0]*x + [1]");
      functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[m],l+2)]->SetParameters(-1/4, 4000);
      functions[Form("fitLine_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[m],l+2)]->SetRange(X_MIN[3*m+l], X_MAX[3*m+l]);
      functions[Form("fitGaus_Alpha_%i_GeV_shStart%i",ENERGIES[m],l+2)] = new TF1("gaus", "gaus");
      functions[Form("fitGaus_Alpha_%i_GeV_shStart%i",ENERGIES[m],l+2)]->SetParameters(1800, 12*0.001, 4*0.001);
      functions[Form("fitGaus_Alpha_%i_GeV_shStart%i",ENERGIES[m],l+2)]->SetRange(A_MIN[3*m+l]*0.001, A_MAX[3*m+l]*0.001);
    }
  }

}

void fillQDCPLots(std::map<std::string, TFile*> &inFiles, std::map<std::string, TH1*> &histos, std::map<std::string, TGraphErrors*> &graphs) {
  for (int i{0}; i<N_ENERGIES; ++i) {
    histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV",ENERGIES[i])] = (TH2D*)((TH2D*)inFiles[Form("Calib_%i_GeV", ENERGIES[i])]->Get("GuilCut_QDCUS_vs_QDCScifi"))->Clone(Form("QDCUS_vs_QDCScifi_Calib_%i_GeV", ENERGIES[i]));
    histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV",ENERGIES[i])] = (TH2D*)((TH2D*)inFiles[Form("Test_%i_GeV", ENERGIES[i])]->Get("GuilCut_QDCUS_vs_QDCScifi"))->Clone(Form("QDCUS_vs_QDCScifi_Test_%i_GeV", ENERGIES[i]));
    for (int j{0}; j<N_SH_STARTS; ++j) {
      histos[Form("Calib_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH2D*)((TH2D*)inFiles[Form("Calib_%i_GeV", ENERGIES[i])]->Get(Form("GuilCut_QDCUS_vs_QDCScifi_ShStart_st%i",j+2)))->Clone(Form("QDCUS_vs_QDCScifi_Calib_%i_GeV_st%i", ENERGIES[i], j+2));
      histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)] = (TH2D*)((TH2D*)inFiles[Form("Test_%i_GeV", ENERGIES[i])]->Get(Form("GuilCut_QDCUS_vs_QDCScifi_ShStart_st%i",j+2)))->Clone(Form("QDCUS_vs_QDCScifi_Test_%i_GeV_st%i", ENERGIES[i], j+2));

      double scifiQDC = ((TH1F*)inFiles[Form("Calib_%i_GeV", ENERGIES[i])]->Get(Form("GuilCut_Shower_SciFi_QDC_shStart%i",j+2)))->GetMean();
      double scifiQDCerr = ((TH1F*)inFiles[Form("Calib_%i_GeV", ENERGIES[i])]->Get(Form("GuilCut_Shower_SciFi_QDC_shStart%i",j+2)))->GetStdDev();

      TH2D* h = (TH2D*)inFiles[Form("Calib_%i_GeV", ENERGIES[i])]->Get(Form("GuilCut_QDCUS_vs_QDCScifi_ShStart_st%i",j+2));
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

double computeK(std::map<std::string, TFile*> &inFiles, std::map<std::string, TH1*> &histos, std::map<std::string, TF1*> &functions, std::map<std::string, TGraphErrors*> &graphs) {
  double sum_k{0};
  for (int i{0}; i<N_ENERGIES; ++i) {
    for (int j{0}; j<N_SH_STARTS; ++j) {     
      // fit previously filled scatter plots
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

void testEnergyResolution(std::map<std::string, TFile*> &inFiles, std::map<std::string, TH1*> &histos, std::map<std::string, TGraphErrors*> &graphs, const double k, const double alpha) {
  for (int i{0}; i<N_ENERGIES; ++i) {
    for (int j{0}; j<N_SH_STARTS; ++j) { 
      int upper_x_bin = histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->GetNbins();
      int upper_y_bin = histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetYaxis()->GetNbins();
      for (int binX = 1; binX <= upper_x_bin; ++binX) {
        for (int binY = 1; binY <= upper_y_bin; ++binY) {
          // Access the bin content
          double binContent = histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetBinContent(binX, binY);
          if (binContent>0) {
            double scifiQDC = histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetYaxis()->GetBinCenter(binY);
            double usQDC = histos[Form("Test_USQDC_vs_SciFiQDC_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetXaxis()->GetBinCenter(binX);
            double recEn = k*scifiQDC + alpha*usQDC;
            histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[i],j+2)]->Fill(recEn,binContent);
            if (recEn>0) {
              histos[Form("Scifi/E_reco vs starting wall %i GeV",ENERGIES[i])]->Fill(j+1, scifiQDC*k/recEn);
            }            
          }
        }
      }
      double Ereco = histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetMean();
      double ErecoErr = histos[Form("Ereco_%i_GeV_shStart%i",ENERGIES[i],j+2)]->GetStdDev();

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

void drawOnCanvases(std::map<std::string, TH1*> &histos, std::map<std::string, TF1*> &functions, std::map<std::string, TMultiGraph*> &multigraphs, std::map<std::string, TCanvas*> &canvases) {
  for (auto const& mg : multigraphs) {
    canvases[mg.first] = new TCanvas(mg.first.c_str(),mg.first.c_str(),1980,1020);
    canvases[mg.first]->cd();
    mg.second->Draw("ape");
    canvases[mg.first]->BuildLegend();
  }
}

void calibration() {
  std::map<std::string, TFile*> inFiles;
  std::map<std::string, TH1*> histos;
  std::map<std::string, TF1*> functions;
  std::map<std::string, TGraphErrors*> graphs;
  std::map<std::string, TMultiGraph*> multigraphs;
  std::map<std::string, TCanvas*> canvases;


  TFile* outputFile = new TFile("calibrations.root","recreate");
  definePlots(histos, functions, graphs, multigraphs);
  openFiles(inFiles);

  fillQDCPLots(inFiles, histos, graphs);
  double k = computeK(inFiles, histos, functions, graphs);
  double alpha = computeAlpha(inFiles, histos, functions, graphs, k);
  std::cout<<"mean K:\t"<<k<<std::endl;
  std::cout<<"mean Alpha:\t"<<alpha<<std::endl;
  testEnergyResolution(inFiles, histos, graphs, k, alpha);
  drawOnCanvases(histos, functions, multigraphs, canvases);

  outputFile->cd();

  for (auto const& h : histos) {
    h.second->Write();
  }
  for (auto const& c : canvases) {
    c.second->Write();
  }
  
  //outputFile->Write();
  outputFile->Close();

}