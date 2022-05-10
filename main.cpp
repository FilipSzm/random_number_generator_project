#include <iostream>
#include <cmath>
#include <algorithm>
#include <asa/asa091.cpp>
#include <cdflib/cdflib.cpp>


using namespace std;


class Generator {
private:
    unsigned int index;
    unsigned int lower_mask, upper_mask;
    unsigned int* MT;

public:
    explicit Generator (unsigned int seed = 5489u) {
        MT = new unsigned int[624];
        index = 625;
        lower_mask = (unsigned long)(1 << 31) - 1;
        upper_mask = (~lower_mask) & 0xffffffffu;
        init (seed);
    }

    void init (unsigned int seed) {
        index = 624;
        MT[0] = seed;

        for (int i = 1; i < 624; i++) {
            MT[i] = (unsigned long long)(1812433253 * (MT[i - 1] ^ (MT[i - 1] >> 30)) + 1);
        }
    }
    void twist () {
        for (int i = 0; i < 624; i++) {
            unsigned int x = (MT[i] & upper_mask) + (MT[(i + 1) % 624] & lower_mask);
            unsigned long xA = x >> 1;
            if ((x % 2) != 0) {
                xA ^= 0x9908b0dfu;
            }
            MT[i] = MT[(i + 397) % 624] ^ xA;
        }
        index = 0;
    }

    unsigned int G (unsigned long left = 0, unsigned long right = 0xffffffffu) {
        if (index >= 624) {
            if (index > 624) init(5489u);
            twist();
        }
        unsigned long long y = MT[index];
        y = y ^ ((y >> 11) & 0xffffffffu);
        y = y ^ ((y >> 7) & 0x9d2c5680u);
        y = y ^ ((y >> 15) & 0xefc60000u);
        y = y ^ (y >> 18);
        index++;
        if (right == 0xffffffffu && left == 0) {
            return y & 0xffffffffull;
        }
        return (unsigned long)(left + (y % (right - left + 1)) & 0xffffffffu);
    }

    double J (double a = 0, double b = 1) {
        unsigned int x = G();
        double w = (double) x / (double)0xffffffff;
        return b * w + a;
    }
    int B (double p = .5) {
        if (J() <= p) return 1;
        return 0;
    }
    unsigned int D (unsigned int n, double p) {
        unsigned int x = 0;
        for (int i = 0; i < n; i++) {
            if (J() <= p) x++;
        }
        return x;
    }
    int P (double lambda) {
        int k = -1;
        long double p = 1, L = exp(-lambda);
        do {
            k++;
            p *= J();
        } while (p > L);
        return k;
    }
    double W (double lambda = 1) {
        return -log(J()) / lambda;
    }
    double N (double mu = 0, double sigma = 1) {
        double x = sqrt(-2 * log(J())) * cos(2 * M_PI * J());
        return x * sigma + mu;
    }

    ~Generator () {
        delete[] MT;
    }
};

void runs_test (Generator& gen, int n, bool* out) {
    auto* arrG = new unsigned int[n];
    auto* arrJ = new double [n];
    auto* arrB = new int[n];
    bool before;
    for (int i = 0; i < n; i++) arrG[i] = gen.G();
    for (int i = 0; i < n; i++) arrJ[i] = gen.J();
    for (int i = 0; i < n; i++) arrB[i] = gen.B();

    auto* Gsort = new unsigned int[n];
    copy(arrG, arrG + n, Gsort);
    sort(Gsort, Gsort + n);
    unsigned int Gmed = Gsort[n / 2];

    auto* Jsort = new double[n];
    copy(arrJ, arrJ + n, Jsort);
    sort(Jsort, Jsort + n);
    double Jmed = Jsort[n / 2];

    int KG = 1;
    before = (arrG[0] >= Gmed);
    for (int i = 1; i < n; i++) {
        bool curr = (arrG[i] >= Gmed);
        if (curr != before) KG++;
        before = curr;
    }
    int KJ = 1;
    before = (arrJ[0] >= Jmed);
    for (int i = 1; i < n; i++) {
        bool curr = (arrJ[i] >= Jmed);
        if (curr != before) KJ++;
        before = curr;
    }
    int KB = 1;
    before = (arrB[0] > 0);
    for (int i = 1; i < n; i++) {
        bool curr = (arrB[i] > 0);
        if (curr != before) KB++;
        before = curr;
    }

    double N1G = 0, N2G = 0;
    for (int i = 0; i < n; i++) {
        if (arrG[i] >= Gmed) N2G += 1;
        else N1G +=1;
    }
    double N1J = 0, N2J = 0;
    for (int i = 0; i < n; i++) {
        if (arrJ[i] >= Jmed) N2J += 1;
        else N1J += 1;
    }
    double N1B = 0, N2B = 0;
    for (int i = 0; i < n; i++) {
        if (arrB[i] > 0) N2B += 1;
        else N1B += 1;
    }
    double dn = n;
    double EG = ((2. * N1G * N2G) / dn) + 1.;
    double DG = (2. * N1G * N2G * (2. * N1G * N2G - dn)) / ((dn - 1.) * dn * dn);
    double ZG = (double)(KG - EG) / sqrt(DG);
    double EJ = ((2. * N1J * N2J) / dn) + 1.;
    double DJ = (2. * N1J * N2J * (2. * N1J * N2J - dn)) / ((dn - 1.) * dn * dn);
    double ZJ = (double)(KJ - EJ) / sqrt(DJ);
    double EB = ((2. * N1B * N2B) / dn) + 1.;
    double DB = (2. * N1B * N2B * (2. * N1B * N2B - dn)) / ((dn - 1.) * dn * dn);
    double ZB = (double)(KB - EB) / sqrt(DB);

    if (fabs(ZG) <= 1.96) out[0] = true;
    else out[0] = false;
    if (fabs(ZJ) <= 1.96) out[1] = true;
    else out[1] = false;
    if (fabs(ZB) <= 1.96) out[2] = true;
    else out[2] = false;
    //cout << ZG << " " << ZJ << " " << ZB << endl;

    delete[] arrG;
    delete[] arrJ;
    delete[] arrB;
    delete[] Gsort;
    delete[] Jsort;
}


double chi_square_ppf(double v, double p) {
    int fx;
    double temp = ppchi2(p, v, lgamma (v / 2.0), &fx);
    return temp;
}

bool Chi_Square_G (Generator& gen, int n) {
    int k = pow(n, 4./5);
    double binMax = (double)0xffffffff / k;
    auto* rand = new double [n];
    for (int i = 0; i < n; i++) rand[i] = gen.G();
    auto* binObs = new int[k];
    auto* binExp = new double[k];
    for (int i = 0; i < k; i++) {
        binObs[i] = 0;
    }
    for (int i = 0; i < n; i++) {
        binObs[(int)(rand[i] / binMax)]++;
    }
    for (int i = 0; i < k; i++) {
        binExp[i] = (double)n / k;
    }
    double chi_squared = 0;
    for (int i = 0; i < k; i++) {
        if (i + 2 == k && binObs[i + 1] < 5) {
            i++;
            binObs[i] += binObs[i - 1];
            binExp[i] += binExp[i - 1];
            k--;
        }
        if (binObs[i] >= 5 || i == k - 1) {
            chi_squared += (binObs[i] - binExp[i]) * (binObs[i] - binExp[i]) / binExp[i];
        }
        else if (i + 2 < k){
            binObs[i + 1] += binObs[i];
            binExp[i + 1] += binExp[i];
            k--;
        }
    }


    delete[] rand;
    delete[] binObs;
    delete[] binExp;
    if (chi_squared <= chi_square_ppf(k - 1, 0.95 /*.975*/) /*&& chi_squared >= chi_square_ppf(k - 1, .025)*/) return true;
    return false;
}

bool Chi_Square_J (Generator& gen, int n) {
    int k = 2 * pow(n, 4./5);
    double binMax = 1. / k;
    auto* rand = new double [n];
    for (int i = 0; i < n; i++) rand[i] = gen.J();
    auto* binObs = new int[k];
    auto* binExp = new double[k];
    for (int i = 0; i < k; i++) {
        binObs[i] = 0;
    }
    for (int i = 0; i < n; i++) {
        binObs[(int)(rand[i] / binMax)]++;
    }
    for (int i = 0; i < k; i++) {
        binExp[i] = (double)n / k;
    }
    long double chi_squared = 0;
    for (int i = 0; i < k; i++) {
        if (i + 2 == k && binObs[i + 1] < 5) {
            i++;
            binObs[i] += binObs[i - 1];
            binExp[i] += binExp[i - 1];
            k--;
        }
        if (binObs[i] >= 5 || i == k - 1) {
            chi_squared += (binObs[i] - binExp[i]) * (binObs[i] - binExp[i]) / binExp[i];
        }
        else if (i + 2 < k){
            binObs[i + 1] += binObs[i];
            binExp[i + 1] += binExp[i];
            k--;
        }
    }

    delete[] rand;
    delete[] binObs;
    delete[] binExp;
    if (chi_squared <= chi_square_ppf(k - 1, 0.95 /*.975*/) /*&& chi_squared >= chi_square_ppf(k - 1, .025)*/) return true;
    return false;
}


double bin_cdf (double n, double  p, double a) {
    double temp = pow(p, a) * pow(1 - p, n - a);
    temp /= ((n + 1) * beta(n - a + 1., a + 1.));
    return temp;
}

bool Chi_Square_D (Generator& gen, int n, double p) {
    int k = sqrt(n);
    auto* rand = new unsigned int[n];
    for (int i = 0; i < n; i++) rand[i] = gen.D(k - 1, p);
    auto* binObs = new int[k];
    auto* binExp = new double[k];
    for (int i = 0; i < k; i++) {
        binObs[i] = 0;
    }
    for (int i = 0; i < n; i++) {
        binObs[rand[i]]++;
    }
    for (int i = 0; i < k; i++) {
        binExp[i] = n * bin_cdf(k - 1, p, i);
    }

    long double chi_squared = 0;
    for (int i = 0; i < k; i++) {
        if (i + 2 == k && binObs[i + 1] < 5) {
            i++;
            binObs[i] += binObs[i - 1];
            binExp[i] += binExp[i - 1];
            k--;
        }
        if (binObs[i] >= 5 || i == k - 1) {
            chi_squared += (binObs[i] - binExp[i]) * (binObs[i] - binExp[i]) / binExp[i];
        }
        else if (i + 2 < k){
            binObs[i + 1] += binObs[i];
            binExp[i + 1] += binExp[i];
            k--;
        }
    }

    delete[] rand;
    delete[] binObs;
    delete[] binExp;
    if (chi_squared <= chi_square_ppf(k - 2, 0.95 /*.975*/) /*&& chi_squared >= chi_square_ppf(k - 2, .025)*/) return true;
    return false;
}

double factorial (int x) {
    if (x == 0) return 1;
    return x * factorial(x - 1);
}

bool Chi_Square_P (Generator& gen, int n, double lambda) {
    int k = (2. * lambda) + 4;

    auto* rand = new int [n];
    for (int i = 0; i < n; i++) rand[i] = gen.P(lambda);
    auto* binObs = new int[k];
    auto* binExp = new double[k];
    for (int i = 0; i < k; i++) {
        binObs[i] = 0;
    }
    int n2 = n;
    for (int i = 0; i < n; i++) {
        if (rand[i] >= k) {
            n2--;
        }
        else {
            binObs[rand[i]]++;
        }
    }
    n = n2;
    for (int i = 0; i < k; i++) {
        binExp[i] = (double)n * (pow(lambda, (double)i) / factorial(i) * exp(-lambda));
    }

    double chi_squared = 0;
    for (int i = 0; i < k; i++) {
        if (i + 2 == k && binObs[i + 1] < 5) {
            i++;
            binObs[i] += binObs[i - 1];
            binExp[i] += binExp[i - 1];
            k--;
        }
        if (binObs[i] >= 5 || i == k - 1) {
            chi_squared += (binObs[i] - binExp[i]) * (binObs[i] - binExp[i]) / binExp[i];
        }
        else if (i + 2 < k){
            binObs[i + 1] += binObs[i];
            binExp[i + 1] += binExp[i];
            k--;
        }
    }

    delete[] rand;
    delete[] binObs;
    delete[] binExp;
    if (chi_squared <= chi_square_ppf(k - 2, 0.95 /*.975*/) /*&& chi_squared >= chi_square_ppf(k - 2, .025)*/) return true;
    return false;
}

bool Chi_Square_W (Generator& gen, int n, double lambda) {
    int k = pow(n, 4./5);

    auto* rand = new double [n];
    for (int i = 0; i < n; i++) rand[i] = gen.W(lambda);
    double min = rand[0], max = rand[0];
    auto* binObs = new int[k];
    auto* binExp = new double[k];
    for (int i = 0; i < n; i++) {
        if (rand[i] < min) min = rand[i];
        if (rand[i] > max) max = rand[i];
    }
    double binMax = (max - min) / k;
    for (int i = 0; i < k; i++) {
        binObs[i] = 0;
    }
    for (int i = 0; i < n; i++) {

        binObs[(int)((rand[i] - min) / binMax)]++;

    }
    for (int i = 0; i < k; i++) {
        binExp[i] = (double)n * (exp(-lambda * (min + binMax * (double)i)) - exp(-lambda * (min + binMax * (double)(i + 1))));
    }

    double chi_squared = 0;
    for (int i = 0; i < k; i++) {
        if (i + 2 == k && binObs[i + 1] < 5) {
            i++;
            binObs[i] += binObs[i - 1];
            binExp[i] += binExp[i - 1];
            k--;
        }
        if (binObs[i] >= 5 || i == k - 1) {
            chi_squared += (binObs[i] - binExp[i]) * (binObs[i] - binExp[i]) / binExp[i];
        }
        else if (i + 2 < k){
            binObs[i + 1] += binObs[i];
            binExp[i + 1] += binExp[i];
            k--;
        }
    }

    delete[] rand;
    delete[] binObs;
    delete[] binExp;
    if (chi_squared <= chi_square_ppf(k - 2, 0.95 /*.975*/) /*&& chi_squared >= chi_square_ppf(k - 2, .025)*/) return true;
    return false;
}

double normalCDF(double value) {
    return 0.5 * erfc(-value * M_SQRT1_2);
}

bool Chi_Square_N (Generator& gen, int n) {
    int k = pow(n, 4./5);

    auto* rand = new double [n];
    for (int i = 0; i < n; i++) rand[i] = gen.N();
    double min = rand[0], max = rand[0];
    auto* binObs = new int[k];
    auto* binExp = new double[k];
    for (int i = 0; i < n; i++) {
        if (rand[i] < min) min = rand[i];
        if (rand[i] > max) max = rand[i];
    }
    double binMax = (max - min) / k;
    for (int i = 0; i < k; i++) {
        binObs[i] = 0;
    }
    for (int i = 0; i < n; i++) {

        binObs[(int)((rand[i] - min) / binMax)]++;

    }
    for (int i = 0; i < k; i++) {
        binExp[i] = (double)n * (normalCDF(min + binMax * (double)(i + 1)) - normalCDF(min + binMax * (double)i));
    }

    double chi_squared = 0;
    for (int i = 0; i < k; i++) {
        if (i + 2 == k && binObs[i + 1] < 5) {
            i++;
            binObs[i] += binObs[i - 1];
            binExp[i] += binExp[i - 1];
            k--;
        }
        if (binObs[i] >= 5 || i == k - 1) {
            chi_squared += (binObs[i] - binExp[i]) * (binObs[i] - binExp[i]) / binExp[i];
        }
        else if (i + 2 < k) {
            binObs[i + 1] += binObs[i];
            binExp[i + 1] += binExp[i];
            k--;
        }
    }

    delete[] rand;
    delete[] binObs;
    delete[] binExp;
    if (chi_squared <= chi_square_ppf(k - 2, 0.95 /*.975*/) /*&& chi_squared >= chi_square_ppf(k - 2, .025)*/) return true;
    return false;
}

bool Chi_Square_B (Generator& gen, int n, double p) {
    int k = 2;

    auto* rand = new int[n];
    for (int i = 0; i < n; i++) rand[i] = gen.B(p);
    auto* binObs = new int[k];
    auto* binExp = new double[k];
    for (int i = 0; i < k; i++) {
        binObs[i] = 0;
    }
    for (int i = 0; i < n; i++) {
        binObs[rand[i]]++;
    }
    binExp[0] = (double)n * (1. - p);
    binExp[1] = (double)n * (p);

    double chi_squared = 0;
    for (int i = 0; i < k; i++) {
        if (i + 2 == k && binObs[i + 1] < 5) {
            i++;
            binObs[i] += binObs[i - 1];
            binExp[i] += binExp[i - 1];
            k--;
        }
        if (binObs[i] >= 5 || i == k - 1) {
            chi_squared += (binObs[i] - binExp[i]) * (binObs[i] - binExp[i]) / binExp[i];
        }
        else if (i + 2 < k) {
            binObs[i + 1] += binObs[i];
            binExp[i + 1] += binExp[i];
            k--;
        }
    }

    delete[] rand;
    delete[] binObs;
    delete[] binExp;
    if (chi_squared <= chi_square_ppf(k - 1, 0.95 /*.975*/) /*&& chi_squared >= chi_square_ppf(k - 2, .025)*/) return true;
    return false;
}

int main() {
    auto* gen = new Generator();
    bool* out = new bool[3];

    runs_test(*gen, 10000, out);
    cout << out[0] << " " << out[1] << " " << out[2] << endl;

    cout << endl;

    cout << Chi_Square_G(*gen, 1000) << " ";
    cout << Chi_Square_J(*gen, 1000) << " ";
    cout << Chi_Square_D(*gen, 10000, 1./2) << " ";
    cout << Chi_Square_P(*gen, 10000, 1./2) << " ";
    cout << Chi_Square_W(*gen, 1000, 1) << " ";
    cout << Chi_Square_N(*gen, 10000) << " ";
    cout << Chi_Square_B(*gen, 1000, 1./5) << " ";

    cout << endl << endl;

    for (int i = 0; i < 10000; i++) {
        cout << gen->N() << endl;
    }
    cout << endl;


    delete[] out;
    delete gen;
    return 0;
}



