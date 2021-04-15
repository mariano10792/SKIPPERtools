#include "TChain.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TH1.h"
void run(){
gROOT->SetBatch(1);
TChain data("pvalues");
data.Add("pvalues_sim_hits_corr_proc_skp_moduleC40_41-ssc16_17-lta20_60_TEMP135K-run*_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR600_1_*0.*t");
TCanvas c("");
c.Divide(2,2);
c.cd(1);
data.Draw("pvalues>>h(10,0,1)","ohdu==1");
TH1F * h = new TH1F();
h = (TH1F*)gDirectory->Get("h");
char firstbin[40];
sprintf(firstbin,  "p(p-value<0.1) = %3.2f",h->GetBinContent(1)/h->GetEntries());
TPaveLabel *px1 = new TPaveLabel(0.62,0.15,0.83,0.24,firstbin,"brNDC");
px1->SetBorderSize(0);px1->SetTextFont(42); px1->SetTextSize(0.5);px1->SetTextAlign(22);px1->SetTextColor(1); 
px1->Draw();
c.cd(2);
data.Draw("pvalues>>h2(10,0,1)","ohdu==2");
TH1F * h2 = new TH1F();
h2 = (TH1F*)gDirectory->Get("h2");
char secondbin[40];
sprintf(secondbin,  "p(p-value<0.1) = %3.2f",h2->GetBinContent(1)/h2->GetEntries());
TPaveLabel *px2 = new TPaveLabel(0.62,0.15,0.83,0.24,secondbin,"brNDC");
px2->SetBorderSize(0);px2->SetTextFont(42); px2->SetTextSize(0.5);px2->SetTextAlign(22);px2->SetTextColor(1); 
px2->Draw();
c.cd(3);
data.Draw("pvalues>>h3(10,0,1)","ohdu==3");
TH1F * h3 = new TH1F();
h3 = (TH1F*)gDirectory->Get("h3");
char thirdbin[40];
sprintf(thirdbin,  "p(p-value<0.1) = %3.2f",h3->GetBinContent(1)/h3->GetEntries());
TPaveLabel *px3 = new TPaveLabel(0.62,0.15,0.83,0.24,thirdbin,"brNDC");
px3->SetBorderSize(0);px3->SetTextFont(42); px3->SetTextSize(0.5);px3->SetTextAlign(22);px3->SetTextColor(1); 
px3->Draw();
c.cd(4);
data.Draw("pvalues>>h4(10,0,1)","ohdu==4");
TH1F * h4 = new TH1F();
h4 = (TH1F*)gDirectory->Get("h4");
char fourthbin[40];
sprintf(fourthbin,  "p(p-value<0.1) = %3.2f",h4->GetBinContent(1)/h4->GetEntries());
TPaveLabel *px4 = new TPaveLabel(0.62,0.15,0.83,0.24,fourthbin,"brNDC");
px4->SetBorderSize(0);px4->SetTextFont(42); px4->SetTextSize(0.5);px4->SetTextAlign(22);px4->SetTextColor(1); 
px4->Draw();


c.SaveAs("Module1halo5run9.pdf");
}
