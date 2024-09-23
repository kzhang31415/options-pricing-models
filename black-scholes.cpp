#include <cmath>
#include <vector>
#include <iostream>
#include <numbers>
#include <random>
using namespace std;
#define PI 3.14159265

//Compile with g++ black-scholes.cpp --std=c++17 -o black-scholes on ARM Mac 
class BlackScholesModel{
    public:
        double S0, K, r, sigma;
        double d1, d2;
        BlackScholesModel(double S0, double K, double r, double sigma){
            this->S0 = S0;
            this->K = K;
            this->r = r;
            this->sigma = sigma;
        }

        void update(double S0, double K, double r, double sigma){
            this->S0 = S0;
            this->K = K;
            this->r = r;
            this->sigma = sigma;
        }

        double getCallPrice(double T){
            double d1 = (log(S0/K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
            double d2 = d1 - sigma * sqrt(T);
            return S0 * normalCDF(d1) - K * exp(-r * T) * normalCDF(d2);
        }
        
        double getPutPrice(double T){
            double d1 = (log(S0/K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
            double d2 = d1 - sigma * sqrt(T);
            return K * exp(-r * T) * normalCDF(-d2) - S0 * normalCDF(-d1);
        }

    private:
        double normalCDF(double x){
            return 0.5 * erfc(-x/sqrt(2));
        }
};

class MonteCarloModel{
    public: 
        double S0, K, r, sigma; int n;
        MonteCarloModel(double S0, double K, double r, double sigma, int n){
            this->S0 = S0;
            this->K = K;
            this->r = r;
            this->sigma = sigma;
            this->n = n;
        }

        void update(double S0, double K, double r, double sigma, int n){
            this->S0 = S0;
            this->K = K;
            this->r = r;
            this->sigma = sigma;
            this->n = n;
        }

        vector<double> priceSimulation(double T){
            vector<double> prices(n+1, 0.0);
            prices[0] = S0;
            for(int t = 0; t < n; t++){
                double z = gaussian();
                double S_t = S0 * exp((r - 0.5 * sigma * sigma) * t + sigma * sqrt(t) * z);
                prices[t+1] = S_t;
            }
            return prices;
        }

        double getCallPrice(double T){
            double sum = 0;
            for(int i = 0; i < n; i++){
                double St = S0*exp((r - 0.5*sigma*sigma)*T + sigma*sqrt(T)*gaussian());
                sum += max(St - K, 0.0);
            }
            return exp(-r * T) * (sum / n);
        }

        double getPutPrice(double T){
            double sum = 0;
            for(int i = 0; i < n; i++){
                double St = S0*exp((r - 0.5*sigma*sigma)*T + sigma*sqrt(T)*gaussian());
                sum += max(K - St, 0.0);
            }
            return exp(-r * T) * (sum / n);
        }

    private:
        double gaussian() {
            // Box-Muller transform 
            double x = 0.0;
            double y = 0.0;
            double euclid_sq = 0.0;
            do {
                x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
                y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
                euclid_sq = x*x + y*y;
            } while (euclid_sq >= 1.0);

            return x*sqrt(-2*log(euclid_sq)/euclid_sq);
        }
};

int main(){
    //Example: S0 = 100, K = 100, r = 0.05, sigma = 0.2, T = 1, n = 10000
    double S0, K, r, sigma, T; int n;
    cin >> S0 >> K >> r >> sigma >> T >> n;
    BlackScholesModel bsModel(S0, K, r, sigma);
    MonteCarloModel mcModel(S0, K, r, sigma, n);
    cout << "Exact Call: " << bsModel.getCallPrice(T) << "\n";
    cout << "Exact Put: " << bsModel.getPutPrice(T) << "\n";
    cout << "Monte Carlo Call: " << mcModel.getCallPrice(T) << "\n";
    cout << "Monte Carlo Put: " << mcModel.getPutPrice(T) << "\n";
    return 0;
}
