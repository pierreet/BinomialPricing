
#include <cmath>
#include <algorithm>
#include <vector>


double callEU(double s, double p, double r, double u, double d, int n){

    double Rinv = 1/(1+r);
    double q = ((1+r)-d)/(u-d);
    vector<double> v(n+1);
    vector<double> call(n+1);

    for(int i=0; i<=n; ++i)
        v[i]=pow(u, i)*pow(d, n-i)*p;

    for(int i=0; i<=n; ++i)
        call[i]=max(0.0, (v[i]-s));

    for(int i=n-1; i>=0; --i){
        for(int j=0; j<=i; ++j)
            call[j] = (q*call[j+1]+(1-q)*call[j])*Rinv;
    }

    return call[0];
}

double putEU(double s, double p, double r, double u, double d, int n){

    double Rinv = 1/(1+r);
    double q = ((1+r)-d)/(u-d);
    vector<double> v(n+1);
    vector<double> put(n+1);

    for(int i=0; i<=n; ++i)
        v[i]=pow(u, i)*pow(d, n-i)*p;

    for(int i=0; i<=n; ++i)
        put[i]=max(0.0, (s-v[i]));

    for(int i=n-1; i>=0; --i){
        for(int j=0; j<=i; ++j)
            put[j] = (q*put[j+1]+(1-q)*put[j])*Rinv;
    }

    resultat = put[0];
}

double putUS(double s, double p, double r, double u, double d, int n){

    double Rinv = 1/(1+r);
    double q = ((1+r)-d)/(u-d);
    vector<double> v(n+1);
    vector<double> put(n+1);

    for(int i=0; i<=n; ++i)
        v[i]=pow(u, i)*pow(d, n-i)*p;

    for(int i=0; i<=n; ++i)
        put[i]=max(0.0, (s-v[i]));

    for(int i=n-1; i>=0; --i){
        for(int j=0; j<=i; ++j)
            put[j] = max((q*put[j+1]+(1-q)*put[j])*Rinv, (s-(pow(u, j)*pow(d, i-j)*p)));

    }

    return put[0];
}

double callUS(double s, double p, double r, double u, double d, int n){

    double Rinv = 1/(1+r);
    double q = ((1+r)-d)/(u-d);
    vector<double> v(n+1);
    vector<double> call(n+1);

    for(int i=0; i<=n; ++i)
        v[i]=pow(u, i)*pow(d, n-i)*p;

    for(int i=0; i<=n; ++i)
        call[i]=max(0.0, (v[i]-s));

    for(int i=n-1; i>=0; --i){
        for(int j=0; j<=i; ++j)
            call[j] = max((q*call[j+1]+(1-q)*call[j])*Rinv, ((pow(u, j)*pow(d, i-j)*p)-s));

    }

    return call[0];
}

