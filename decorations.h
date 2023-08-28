/*void decorate_canvas(TCanvas *c1){
  TPad *mainPad;
  
  }*/

TPad *create_mainpad(){
  TPad *pad = new TPad("", "", 0.0, 0.25, 1.0, 1.0);
  pad->SetLeftMargin(0.1);
  pad->SetRightMargin(0.25);
  pad->SetTopMargin(0.09);
  pad->SetBottomMargin(0.12);
  pad->SetTickx(1);
  pad->SetTicky(1);
  return pad;
}

TPad *create_ratiopad(){
  TPad *pad = new TPad("", "", 0.0, 0.01, 1.0, 0.25);
  pad->SetLeftMargin(0.1);
  pad->SetRightMargin(0.25);
  //pad->SetTopMargin(0.09);
  pad->SetBottomMargin(0.12);
  pad->SetTickx(1);
  pad->SetTicky(1);
  return pad;
}

TLegend *create_legend(){
  TLegend *lg1 = new TLegend(0.9, 0.45, 0.77,0.86);
  lg1->SetTextFont(62);		
  lg1->SetFillStyle(0);
  lg1->SetBorderSize(0);
  lg1->SetTextSize(0.03);
  return lg1;
}

void SetHistoStyle(TH1F *h, TString plotname, int color){
  h->SetLineColor(color);
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.2);
  h->SetTitle("");

  //Y-axis
  h->GetYaxis()->SetTitle("Events");
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(18);
  h->GetYaxis()->SetTitleOffset(1.7);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(15);
  h->GetYaxis()->SetLabelOffset(0.02);
  h->GetYaxis()->CenterTitle();
  //h->GetYaxis()->SetTickLength(-0.02);

  //X-axis
  h->GetXaxis()->SetTitle(plotname);
  h->GetXaxis()->SetTitleFont(43);
  h->GetXaxis()->SetTitleSize(16);
  h->GetXaxis()->SetTitleOffset(1.5);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(15);
  h->GetXaxis()->SetLabelOffset(0.02);
  //h->GetXaxis()->CenterTitle();
  //h->GetXaxis()->SetTickLength(-0.02);
}

void SetHistoStyle_bold(TH1F *h, TString plotname, int color){
  SetHistoStyle(h, plotname,color);
  h->SetLineWidth(2);
}

void PutText(TString jobname){
  float x = 0.1;
  float y = 0.95;
  TString job = "jobname = "+jobname;
  TLatex *t = new TLatex();
  t->SetTextFont(42);
  t->SetTextSize(0.037);
  t->DrawLatex(x, y, "CMS");
  t->SetTextSize(0.03);
  //t->DrawLatex(x + 0.4, y, job);
  t->SetTextFont(72);
  t->SetTextSize(0.04);
  t->DrawLatex(x + 0.07, y, "work in progress");
}
