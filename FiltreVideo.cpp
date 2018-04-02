#include <opencv2/opencv.hpp>
#include "IIRFilter.hpp"

using namespace std;
using namespace cv;

struct ParamVideo{
    int fLow;
    int fHigh;
    int ordre; // Ordre du filtre
    int typeFiltre; //0=lp, 1 bp, 2=hp
    IIRFilter *filter;
    vector<Mat> x;
    vector<Mat> y;
    String nomfenetre;
};


void AjouteGlissiere(String nomGlissiere, String nomFenetre, int minGlissiere, int maxGlissiere, int valeurDefaut, int *valGlissiere, void(*f)(int, void *), void *r = NULL);
void MAJFiltre(int x, void * r);
Mat IIRtemporalFilter(IIRFilter &f, Mat phase, vector<Mat> &ri);

Mat IIRtemporalFilter(IIRFilter &f, Mat phase, vector<Mat> &ri)
{
    Mat tmp;
    tmp = f.b[0] * phase + (ri[0]);
    for (int i = 0; i<ri.size() - 2; i++)
        ri[i] = f.b[i + 1] * phase + (ri[i + 1]) + f.a[i + 1] * tmp;
    ri[ri.size() - 2] = f.b[ri.size() - 1] * phase + f.a[ri.size() - 1] * tmp;
    return tmp;
}

int main(int argc, char *argv[])
{
    VideoCapture v(0);
    if (!v.isOpened())
    {
        cout << "Cannot open video\n";
        return 0;
    }
    ParamVideo p;
    p.nomfenetre = "video";
    namedWindow(p.nomfenetre);
    p.fLow = 0;
    p.fHigh = 20;
    p.typeFiltre = 0;
    p.ordre = 2;
    p.filter = NULL;
    AjouteGlissiere("fLow", p.nomfenetre, 0, 50, p.fLow, &p.fLow, MAJFiltre, &p);
    AjouteGlissiere("fHigh", p.nomfenetre, 0, 50, p.fHigh, &p.fHigh, MAJFiltre, &p);
    AjouteGlissiere("Order", p.nomfenetre, 1, 3, p.ordre, &p.ordre, MAJFiltre, &p);
    vector<double> f = { p.fLow / 100.0,p.fHigh / 100.0 };
    p.filter = new IIRFilter("butt", p.ordre, 1, f);
    int  code = 0;
    Mat frameUSB;
    Mat frame;
    v >> frameUSB;
    frameUSB.convertTo(frame, CV_32F);
    for (int i = 0; i <= 20; i++)
    {
        p.x.push_back(frame.clone());
        p.y.push_back(frame.clone());
    }
    while (code != 27)
    {
        v >> frameUSB;
        frameUSB.convertTo(frame, CV_32F);
        p.x[0] = frame.clone();
        Mat r;
        r = p.filter->b[0] * p.x[0];
        for (int i = 1; i < p.filter->b.size(); i++)
        {
            r += p.filter->b[i] * p.x[i]; 
        }
        for (int i = 1; i < p.filter->a.size(); i++)
        {
            r += p.filter->a[i] * p.y[i];
        }
        for (int i = p.filter->b.size() - 1; i > 0; i--)
        {
            p.x[i] = p.x[i - 1];
        }
        for (int i = p.filter->a.size() - 1; i > 1; i--)
        {
            p.y[i] = p.y[i - 1];
        }
        p.y[1] = r;
        imshow(p.nomfenetre, r/256);
        code = waitKey(1);
    }
}

void MAJFiltre(int x, void * r)
{
    ParamVideo *p= (ParamVideo*) r;
    delete p->filter;
    vector<double> f = { p->fLow / 100.0,p->fHigh / 100.0 };
    p->filter = new IIRFilter("butt", p->ordre, 1, f);
}

void AjouteGlissiere(String nomGlissiere, String nomFenetre, int minGlissiere, int maxGlissiere, int valeurDefaut, int *valGlissiere, void(*f)(int, void *), void *r)
{
    createTrackbar(nomGlissiere, nomFenetre, valGlissiere, 1, f, r);
    setTrackbarMin(nomGlissiere, nomFenetre, minGlissiere);
    setTrackbarMax(nomGlissiere, nomFenetre, maxGlissiere);
    setTrackbarPos(nomGlissiere, nomFenetre, valeurDefaut);
}
