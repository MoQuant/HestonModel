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

double dWT(){
    int num = 10;
    double dw = (rand() % (2*num + 1)) - num;
    return dw/100.0;
}

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

void CYCLONE(ptree df, std::vector<double> & prices, bool results){
    for(ptree::const_iterator it = df.begin(); it != df.end(); ++it){
        if(results == true){
            for(ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt){
                if(jt->first == "c"){
                    prices.push_back(jt->second.get_value<double>());
                }
            }
        }
        if(it->first == "results"){
            CYCLONE(it->second, std::ref(prices), true);
        }
    }
}

ptree JSON(std::string message){
    std::stringstream fp(message);
    ptree result;
    read_json(fp, result);
    return result;
}

std::vector<double> Computer(std::vector<double> close, double dt){
    auto mean = [](std::vector<double> x){
        double total = 0;
        for(auto & t : x){
            total += t;
        }
        return total / (double) x.size();
    };

    auto variance = [&](std::vector<double> x){
        double total = 0;
        double mu = mean(x);
        for(auto & t : x){
            total += pow(t - mu, 2);
        }
        return total / ((double) x.size() - 1);
    };
    
    std::vector<double> result, ror, sim_vol;
    for(int i = 1; i < close.size(); ++i){
        ror.push_back(close[i]/close[i-1] - 1.0);
    }

    int window = 75;
    for(int i = window; i < ror.size(); ++i){
        std::vector<double> hold = {ror.begin() + (i - window), ror.begin() + i};
        double vol = variance(hold);
        sim_vol.push_back(vol);
    }

    double theta = sqrt(mean(sim_vol));
    double sigma = sqrt(variance(sim_vol));

    double top = 0, bot = 0, mu_vol = mean(sim_vol);
    for(int i = 1; i < sim_vol.size(); ++i){
        top += (sim_vol[i] - mu_vol)*(sim_vol[i-1] - mu_vol);
        bot += pow(sim_vol[i] - mu_vol, 2);
    }

    double kappa = -log(top/bot)/dt;

    double drift = mean(ror);

    result = {kappa, theta, sigma, close[close.size() - 1], drift};
    return result;
}

int main()
{
    srand(time(NULL));

    std::string ticker = "MSFT";
    std::string response = GET(ticker);
    std::vector<double> prices;

    CYCLONE(JSON(response), std::ref(prices), false);

    double t = 30;
    int n = 1000;
    int sims = 100;
    double dt = t / (double) n;

    std::vector<double> params = Computer(prices, dt);
    
    int Long = 0, Short = 0;

    for(int k = 0; k < 100; ++k){
        double forecasted = 0;
        for(int i = 0; i < sims; ++i){
            double S0 = params[3];
            double v0 = 0.1;
            for(int j = 0; j < n; ++j){
                double DWT = dWT();
                v0 += params[0]*(params[1] - v0)*dt + params[2]*sqrt(v0)*-0.5*DWT;
                S0 += params[4]*S0*dt + sqrt(v0)*S0*DWT;
            }
            forecasted += S0;
        }

        forecasted /= (double) sims;

        if(params[3] <= forecasted){
            Long += 1;
        } else {
            Short += 1;
        }

        std::cout << "Long: " << Long << "\tShort: " << Short << std::endl;
    }
    
    std::cout << "Todays Stock Price: " << params[3] << std::endl;

    return 0;
}