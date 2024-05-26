#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <curl/curl.h>
#include <math.h>
#include <cmath>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <algorithm>
#include <time.h>

using namespace boost::property_tree;


// Enter polygon key into this function
std::string URL(std::string ticker)
{
    std::string key = "";
    return "https://api.polygon.io/v2/aggs/ticker/" + ticker + "/range/1/day/2021-01-22/2024-05-26?adjusted=true&sort=asc&limit=1000&apiKey=" + key;
}

size_t WriteCallback(void* contents, size_t size, size_t nmemb, void* userp) {
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

std::string GET(std::string ticker) {

    std::string url = URL(ticker);

    CURL* curl;
    CURLcode res;
    std::string readBuffer;

    curl_global_init(CURL_GLOBAL_DEFAULT);
    curl = curl_easy_init();
    if(curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);

        // Perform the request, res will get the return code
        res = curl_easy_perform(curl);
        
        // Check for errors
        if(res != CURLE_OK) {
            std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << std::endl;
        }

        // Always cleanup
        curl_easy_cleanup(curl);
    }

    curl_global_cleanup();
    return readBuffer;
}

void TYPHOON(ptree df, std::vector<double> & prices, bool results){
    for(ptree::const_iterator it = df.begin(); it != df.end(); ++it){
        if(results == true){
            for(ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt){
                if(jt->first == "c"){
                    prices.push_back(jt->second.get_value<double>());
                }
            }
        }
        if(it->first == "results"){
            TYPHOON(it->second, std::ref(prices), true);
        }
    }
}

ptree JSON(std::string resp){
    std::stringstream ss(resp);
    ptree result;
    read_json(ss, result);
    return result;
}

std::vector<double> Computer(std::vector<double> prices, double dt){
    
    std::vector<double> result, sim_vol, ror;

    auto mean = [](std::vector<double> x){
        double total = 0;
        for(auto & i : x){
            total += i;
        }
        return total / (double) x.size();
    };

    auto variance = [&](std::vector<double> x){
        double total = 0;
        double mu = mean(x);
        for(auto & i : x){
            total += pow(i - mu, 2);
        }
        return total / ((double) x.size() - 1);
    };

    for(int i = 1; i < prices.size(); ++i){
        ror.push_back(prices[i]/prices[i-1] - 1.0);
    }

    double price = prices[prices.size() - 1];
    double drift = mean(ror);

    int window = 50;
    for(int i = window; i < ror.size(); ++i){
        std::vector<double> hold = {ror.begin() + (i - window), ror.begin() + i};
        double vol = variance(hold);
        sim_vol.push_back(vol);
    }

    double theta = sqrt(mean(sim_vol));
    double sigma = sqrt(variance(sim_vol));
    double mu_vol = mean(sim_vol);

    double top = 0, bot = 0;
    for(int i = 1; i < sim_vol.size(); ++i){
        top += (sim_vol[i] - mu_vol)*(sim_vol[i-1] - mu_vol);
        bot += pow(sim_vol[i] - mu_vol, 2);
    }

    double kappa = -log(top/bot)/dt;

    result = {kappa, theta, sigma, price, drift};

    return result;
}

double dWT(){
    int num = 10;
    double dw = (rand() % (2*num + 1)) - num;
    return dw/100.0;
}

int main()
{
    srand(time(NULL));

    std::string ticker = "AMZN";
    std::string response = GET(ticker);

    std::vector<double> prices;
    TYPHOON(JSON(response), std::ref(prices), false);

    double t = 30;
    int n = 1500;
    int p = 150;
    double dt = t / (double) n;

    std::vector<double> params = Computer(prices, dt);

    int Go_Long = 0, Go_Short = 0;

    for(int k = 0; k < 100; ++k){
        double forecast_price = 0;
        for(int i = 0; i < p; ++i){
            double S0 = params[3];
            double v0 = 0.3;
            for(int j = 0; j < n; ++j){
                double dw = dWT();
                v0 += params[0]*(params[1] - v0)*dt + params[2]*sqrt(v0)*dw;
                S0 += params[4]*S0*dt + sqrt(v0)*S0*dw;
            }
            forecast_price += S0;
        }
        forecast_price /= (double) p;
        if(params[3] <= forecast_price){
            Go_Long += 1;
        } else {
            Go_Short += 1;
        }
    }

    std::cout << "Long: " << Go_Long << "\tShort: " << Go_Short << "\tCurrent Stock Price: " << params[3] << std::endl;
    

    return 0;
}
