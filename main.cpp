#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>

using namespace std;
//Rastrigin, De jong 5.12, Griewangk 600 ; Rosenbrock 2.048, SHCB 3,2
double infL = -5.12, supL = 5.12;
double infL2 = -2, supL2 =2;
int popSize=100;
int nrGen=1000;
double Precizie = 3;
double pm = 0.01;
double pcx = 0.4;
int numarP = 30;

vector<double, allocator<double>> toD(const bool crom[]);
double Rastrigin(vector<double> prm);
double Griewangk(vector<double> prm);
double Rosenbrock(vector<double> prm);
double Dejong(vector<double> prm);
double SHCB(double prm[]);
int decimal(double prm[]);
double fRand(double min, double max);

double N = (supL - infL) * pow(10, Precizie);
int d = ceil(log2(N));
int sizeCr = d * numarP;


int main()
{
        bool pop[100][sizeCr];
        bool newGen[100][sizeCr];
        double prob[100];
        int ind[100];
        int i, j, ii, k, nrcx;
        int poz;
        double f[100], fit[100], fs[100], maxf, minf;
        srand(time(nullptr));
        bool sortat;
        double aux;
        int aux2;
        bool aux3;
        double minfG = 10000;
        //generare populatie
        for (i = 0; i < 100; i++)
            for (j = 0; j <= sizeCr; j++) {
                pop[i][j] = rand() % 2;
            }

        //Algoritmul genetic-----------------------------------------------------------------------------------
        int crGen = 1;//counter generatii
        while (crGen < nrGen) {
            //cout<<crGen<<". ";
            //probabilitati
            for (i = 0; i < 100; i++)
                ind[i] = i;
            for (i = 0; i < 100; i++)
                prob[i] = fRand(0, 1);

            //Mutation
            for (i = 0; i < 100; i++)
                for (j = 0; j < sizeCr; j++)
                    if (fRand(0, 1) < pm) {
                        pop[i][j] = 1 - pop[i][j];
                    }

            //Crossover
            //C1.sortare prob si reordonare cromozomi
            // bubblesort
            do {
                sortat = true;
                for (i = 0; i < popSize - 1; i++) {
                    if (prob[i] > prob[i + 1]) {
                        aux = prob[i];
                        prob[i] = prob[i + 1];
                        prob[i + 1] = aux;
                        aux2 = ind[i];
                        ind[i] = ind[i + 1];
                        ind[i + 1] = aux2;
                        sortat = false;
                    }
                }
            } while (!sortat);

            //C2.schimb informatie
            for (i = 0; i < popSize; i = i + 2) {
                if ((prob[i] <= pcx) && (prob[i + 1] <= pcx)) {  //daca prob ambilor cromozomi sunt ca pcx, fac schimbul
                    // de la poz generata random pana la final
                    poz = fRand(1, sizeCr - 2);
                    for (ii = poz; ii < sizeCr; ii++) {
                        aux3 = pop[ind[i]][ii];
                        pop[ind[i]][ii] = pop[ind[i + 1]][ii];
                        pop[ind[i + 1]][ii] = aux3;
                    }
                } else break;
            }

            //Selection-prep
            maxf = -10000;
            minf = 10000;//de revazut la fiecare functie cercetata
            for (i = 0; i < popSize; i++) {
                f[i] = Dejong(toD(pop[i]));//-----------------------------------------------eval fc
                if (maxf < f[i]) maxf = f[i];
                if (minf > f[i]) minf = f[i];
            }
            for (i = 0; i < popSize; i++) {
                fit[i] = 1.1 * maxf - f[i];
            }

            //Roata norocului
            fs[0] = fit[0];
            for (i = 1; i < popSize; i++) {
                fs[i] = fs[i - 1] + fit[i];
            }

            //Selection
            for (i = 0; i < popSize; i++) {
                poz = fRand(0, fs[popSize - 1]);
                for (j = 0; j < popSize; j++) {
                    if (poz <= fs[j]) {//aleg cromozomul de pe poz j din pop in NewGen pe poz i
                        //newGen=pop
                        for (k = 0; k < sizeCr; k++)
                            newGen[i][k] = pop[j][k];
                        break;
                    }
                }
            }
            //pop=newGen;
            for (i = 0; i < popSize; i++) {   //newGen->pop
                for (k = 0; k < sizeCr; k++)
                    pop[i][k] = newGen[i][k];
            }
            crGen++;
            if (minfG > minf) minfG = minf;
        }
        cout << minfG << '\n';

    return 0;
}

double decimal(const bool T[])
{
    int i;
    double x = 0;
    for (i = 0; i < d; ++i)
    {
        x = x * 2;
        x = x + T[i];
    }
    x = x / (pow(2, d) - 1);
    x = x * (supL - infL);
    x = x + infL;
    return x;
}

vector<double> toD(const bool crom[])
{
    //va returna un vector de double
    //argumentul functiilor -> poate avea 5,10,15 dim double codificate
    //decodifica sirul de bool
    int i, j, poz;
    bool temp[d];
    double sol[numarP];
    for (i = 0; i < numarP; i++)
    {
        poz = 0;
        for (j = i * d; j < ((i + 1) * d); j++)
        {
            temp[poz] = crom[j];
            poz++;
        }
        sol[i] = decimal(temp);
    }
    int n = sizeof(sol) / sizeof(sol[0]);
    vector<double> dest(sol, sol + n);
    return dest;
}

double fRand(double min, double max)
{
    double f = (double) rand() / RAND_MAX;
    return min + f * (max - min);
}

double Rastrigin(vector<double> prm)
{
    double S = 0;
    int k;
    for (k = 0; k < numarP; k++)
        S = S + prm[k] * prm[k] - 10 * cos(2 * M_PI * prm[k]);
    S=S+numarP*10;
    return S;
}

double Griewangk(vector<double> prm)
{
    double S = 0, P = 1;
    int k;
    for (k = 0; k < numarP; k++)
    {
        S = S + prm[k] * prm[k] / 4000;
        P = P * cos(prm[k] / sqrt(k+1));
    }

    return S-P+1;
}

double shcb(vector<double> prm)
{
    double S=0;
    S=(4-2.1*prm[0]*prm[0]+prm[0]*prm[0]*prm[0]*prm[0]/3)*prm[0]*prm[0]+prm[0]*prm[1]+(-4+4*prm[1]*prm[1])*prm[1]*prm[1];
    return S;
}
double Dejong(vector<double> prm)
{
    double S=0;
    int k;
    for(k=0; k<numarP; k++)
        S=S+prm[k]*prm[k];
    return S;
}

double Rosenbrock(vector<double> prm)
{
    double S = 0;
    int k;
    for (k = 0; k < numarP - 1; k++)
        S = S + 100 * (prm[k + 1] - prm[k] * prm[k]) * (prm[k + 1] - prm[k] * prm[k]) + (1 - prm[k]) * (1 - prm[k]);
    return S;
}

