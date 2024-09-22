#include <cmath>
#include <vector>
#include <iostream>
#include <numbers>
using namespace std;
#define PI 3.14159265

//Compile with g++ black-scholes.cpp --std=c++17 -o black-scholes on ARM Mac 
class BlackScholesModel{
    public:
        double S0, K, r, sigma, T;
        double d1, d2;
        BlackScholesModel(double S0, double K, double r, double sigma, double T){
            this->S0 = S0;
            this->K = K;
            this->r = r;
            this->sigma = sigma;
            this->T = T;
            this->d1 = (log(S0/K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
            this->d2 = d1 - sigma * sqrt(T);
        }

        void update(double S0, double K, double r, double sigma, double T){
            this->S0 = S0;
            this->K = K;
            this->r = r;
            this->sigma = sigma;
            this->T = T;
            this->d1 = (log(S0/K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
            this->d2 = d1 - sigma * sqrt(T);
        }

        double getCallPrice(){
            return S0 * normalCDF(d1) - K * exp(-r * T) * normalCDF(d2);
        }
        
        double getPutPrice(){
            return K * exp(-r * T) * normalCDF(-d2) - S0 * normalCDF(-d1);
        }

    private:
        double normalCDF(double x){
            return 0.5 * erfc(-x/sqrt(2));
        }
};

int main(){
    //Example: S0 = 100, K = 100, r = 0.05, sigma = 0.2, T = 1
    double S0, K, r, sigma, T;
    cin >> S0 >> K >> r >> sigma >> T;
    BlackScholesModel model(S0, K, r, sigma, T);
    cout << "Call price: " << model.getCallPrice() << "\n";
    cout << "Put price: " << model.getPutPrice() << "\n";
    return 0;
}

