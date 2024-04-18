#include <iostream>
#include <vector>
#include "Fit/Fitter.h"


struct Poinnt {
  double x;
  double y;
};

struct RccCanvas {
  int n;
  TCanvas* c;
  std::vector<TPad*>p;
};



RccCanvas SetupCanvas1(TCanvas *canvas) {
  // Create a TCanvas
  //    TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 600);
  RccCanvas rc;
  rc.c=canvas;
  rc.n=0;
  

  // Divide the canvas into two sections
  // canvas->Divide(1, 2);

  // Get the top pad (upper section)
  // canvas->cd(1);

  // Divide the upper section into a 1x1 pad
  TPad *upperPad = new TPad("upperPad", "Upper Pad", 0, 2.0/3.0, 1, 1);
    
  upperPad->Draw();
  rc.p.push_back((TPad*)upperPad); rc.n++;

  // Get the bottom pad (lower section)
  //canvas->cd(2);

  // Divide the lower section into a 2x2 set of pads
  TPad *lowerPad = new TPad("lowerPad", "Lower Pad", 0, 0, 1, 2.0/3.0);
  lowerPad->Divide(2, 2);
  lowerPad->Draw();
  for (int i=0;i<4;i++){
    rc.p.push_back((TPad*)lowerPad->cd(i+1)); rc.n++;
  }
  return rc;
}



double computeEnclosedArea(const std::vector<Poinnt>& points) {
  //the so-called shoelace algorithm.
    int n = points.size();
    double area = 0.0;

    for (int i = 0; i < n; i++) {
      int j = (i + 1) % n; //the n'th item is also the 0th, for a closed loop as needed.
        area += (points[i].x * points[j].y) - (points[j].x * points[i].y);
    }

    area /= 2.0;
    return std::abs(area);
}

Poinnt computeGeometricCenter(const std::vector<Poinnt>& points) {
  //shoelace-like, taken from https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
    int n = points.size();
    double xSum = 0.0;
    double ySum = 0.0;

    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
	//note that the shoelaceAreaTerm is twice the actual area.
	//we catch the missing factor of 2 outside of the loop.
        double shoelaceAreaTerm =( (points[i].x * points[j].y) - (points[j].x * points[i].y));
        xSum += (points[i].x + points[j].x) * shoelaceAreaTerm;
        ySum += (points[i].y + points[j].y) * shoelaceAreaTerm;
    }

    double area = computeEnclosedArea(points);
    Poinnt center;
    center.x = xSum / (6.0 * area);
    center.y = ySum / (6.0 * area);

    return center;
}

void guessCircleParameters(const std::vector<Poinnt>& points, float &radius, float &x, float &y) {
  //getCenter:
  //we could do this by taking the geometric center, but we know we have only a small arc, so that term will be bad.
  int mid=points.size()/2;
  int last=points.size()-1;
  
  // Calculate midpoints
  double mx1 = (points[0].x + points[mid].x) / 2;
  double my1 = (points[0].y + points[mid].y) / 2;
  double mx2 = (points[mid].x + points[last].x) / 2;
  double my2 = (points[mid].y + points[last].y) / 2;

  // Calculate slopes
  double m1 = (points[mid].y - points[0].y) / (points[mid].x - points[0].x);
  double m2 = (points[last].y - points[mid].y) / (points[last].x - points[mid].x);
  //y=m(x-x0)+y0 ==> m=(y-y0)/(x-x0)

  // Calculate perpendicular bisector slopes
  double m1_perp = -1 / m1;
  double m2_perp = -1 / m2;

  // Calculate intersection point
  x = (my2 - my1 + m1_perp * mx1 - m2_perp * mx2) / (m1_perp - m2_perp);
  y = my1 + m1_perp * (x - mx1);
  // m1(x-x1)+y1=m2(x-x2)+y2
  //m1x-m1x1+y1=m2x-m2x2+y2
  //(m1-m2)x=m1x1-m2x2-y1+y2
  //x=(y2-y1+m1x1-m2x2)

  // Calculate radius
  radius = TMath::Sqrt((points[0].x - x) * (points[0].x - x) + (points[0].y - y) * (points[0].y - y));
   printf("p0=(%f,%f)\np1=(%f,%f)\np2=(%f,%f)\n",points[0].x,points[0].y,points[mid].x,points[mid].y,points[last].x,points[last].y);
   printf("mx1=%f,my1=%f,mx2=%f,my2=%f,m1=%f,m2=%f,perp1=%f,perp2=%f,xc=%f,yc=%f,rad=%f\n",
  	 mx1,my1,mx2,my2,m1,m2,m1_perp,m2_perp,x,y,radius);
   return;
}
ROOT::Fit::FitResult iterateCircleCenter(const std::vector<Poinnt>& points) {
  //taken from https://link.springer.com/article/10.1007/s10851-010-0249-8
  //getCenter:
  //we could do this by taking the geometric center, but we know we have only a small arc, so that term will be bad.
  float radius,xc,yc;
  guessCircleParameters(points, radius,xc,yc);

   auto chi2Function = [&](const Double_t *par) {
      //minimisation function computing the sum of squares of residuals
     Int_t np = points.size();
      Double_t f = 0;
       for (Int_t i=0;i<np;i++) {
         Double_t u = points[i].x - par[0];
         Double_t v = points[i].y - par[1];
         Double_t dr = par[2] - std::sqrt(u*u+v*v);
         f += dr*dr;
      }
     //  printf("points=%d, center=(%.1f,%.1f) rad=%.1f, res=%f\n",np,par[0],par[1],par[2],f);
     return f;
   };
   // wrap chi2 function in a function object for the fit
   // 3 is the number of fit parameters (size of array par)
   ROOT::Math::Functor fcn(chi2Function,3);
   ROOT::Fit::Fitter  fitter;
   //printf("starting with xc=%f,yc=%f\n",xc,yc);
   double pStart[3] = {xc+3,yc-3,radius/2};//xc,yc,radius};
   fitter.SetFCN(fcn, pStart);
   fitter.Config().SetMinimizer("Minuit2", "Migrad");
   //   fitter.somethingToSetTheMinimizer
   //fitter.Config().NPar();
   fitter.Config().ParSettings(0).Set("x0",xc,0.02,xc-3,xc+3);
   //fitter.Config().ParSettings(0).SetName("x0");
   //  fitter.Config().ParSettings(0).SetLowerLimit(xc+5);
   // fitter.Config().ParSettings(0).SetUpperLimit(xc+10);
   fitter.Config().ParSettings(1).Set("y0",yc,0.02,yc-3,yc+3);
  //fitter.Config().ParSettings(1).SetName("y0");
   //fitter.Config().ParSettings(1).SetLowerLimit(yc-1);
   //fitter.Config().ParSettings(1).SetUpperLimit(yc+1);
   // fitter.Config().ParSettings(2).Set("R",radius/2);//,0.02,2,6);
   fitter.Config().ParSettings(2).SetName("R");
   //fitter.Config().ParSettings(2).SetUpperLimit(9);
   //printf("about to do the fit...\n");
   //printf("fitter is:%s\n",fitter.Config().MinimizerName().c_str());
   //printf("fitter has %d params\n",fitter.Config().NPar());
   // printf("minimizer has %d params\n",fitter.GetMinimizer()->NFree());
   // do the fit 
   bool ok = fitter.FitFCN();
   if (!ok) {
     Error("line3Dfit","Line3D Fit failed");
   }
   const ROOT::Fit::FitResult & result = fitter.Result();
   result.Print(std::cout);

   return fitter.Result();

}

void opc_to_tree(){
  TFile *input=TFile::Open("petal55.ogp_tree.root");
  TTree *ogpTree=(TTree*)(input->Get("ogpTree"));

  //draw everything
  TCanvas *canvas=new TCanvas(Form("ogp%d",55),Form("ogp%d",55),800,600);
  RccCanvas c1=SetupCanvas1(canvas);
  c1.p[0]->cd();
  ogpTree->Draw("p.x():p.y()");//,"ogpIndex==7","*");

  //cycle over the 'isCircle' sets, which are 7,8,9,10
  
  for (int i=7;i<11;i++){
    c1.p[i-6]->cd();
    int nPoints=ogpTree->Draw("p.x():p.y()",Form("ogpIndex==%d",i),"*");    
    std::vector<Poinnt> points;
    for (int j=0;j<nPoints;j++){
      Poinnt temp;
      temp.x=ogpTree->GetV1()[j];
      temp.y=ogpTree->GetV2()[j];
      points.push_back(temp);
      //points[j].x=temp.x;
      //points[j].y=temp.y;
      
    }
    printf("completed fit, drawing:\n");
    auto result=iterateCircleCenter(points);
    TArc *arc = new TArc(result.Parameter(1),result.Parameter(0),result.Parameter(2));
    arc->SetLineColor(kRed);
    arc->SetFillStyle(0);
    arc->Draw();
    if(0){//block to sanity-check the fit
      float x,y,r;
      guessCircleParameters(points,r,x,y);
      printf("x:%f\ny:%f\nr:%f\n",x,y,r);
      arc = new TArc(y,x,r);//result.Parameter(0),5);
      arc->SetLineColor(kBlue);
      arc->SetFillStyle(0);
      arc->Draw();
      arc = new TArc(points[0].y,points[0].x,0.2);
      arc->SetLineColor(kBlue);
      arc->SetFillStyle(0);
      arc->Draw();
      arc = new TArc(points[points.size()/2].y,points[points.size()/2].x,0.2);
      arc->SetLineColor(kBlue);
      arc->SetFillStyle(0);
      arc->Draw();
      arc = new TArc(points[points.size()-1].y,points[points.size()-1].x,0.2);
      arc->SetLineColor(kBlue);
      arc->SetFillStyle(0);
      arc->Draw();
    }
  }
  c1.c->Draw();
}
