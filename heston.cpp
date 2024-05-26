#include <iostream>
#include <curl/curl.h>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <algorithm>
#include <functional>
#include <time.h>

using namespace boost::property_tree;

std::string key()
{
    return "_qfnK4nx068Hwzfer53IJWlZ_YpF97cM";
}

static size_t CallBack(void *contents, size_t size, size_t nmemb, std::string *output)
{
    size_t totalSize = size * nmemb;
    output->append((char*)contents, totalSize);
    return totalSize;
}

std::string GET(std::string ticker)
{
    std::string url = "https://api.polygon.io/v2/aggs/ticker/" + ticker + "/range/1/day/2021-01-25/2024-05-25?adjusted=true&sort=asc&limit=1000&apiKey=" + key();
    std::string response;
    curl_global_init(CURL_GLOBAL_DEFAULT);
    CURL *curl = curl_easy_init();
    struct curl_slist* headers = NULL;
    headers = curl_slist_append(headers, "Content-Type: application/json");
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, CallBack);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    curl_global_cleanup();
    return response;
}

void Typhoon(ptree data, std::vector<double> & prices)
{
    for(ptree::const_iterator it = data.begin(); it != data.end(); ++it){
        for(ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt){
            if(jt->first == "c"){
                prices.push_back(jt->second.get_value<double>());
            }
        }
    }
}

void Cyclone(ptree data, std::vector<double> & prices)
{
    for(ptree::const_iterator it = data.begin(); it != data.end(); ++it){

        if(it->first == "results"){
            Typhoon(it->second, std::ref(prices));
        }
    }
}

ptree JSON(std::string text)
{
    std::stringstream ss(text);
    ptree result;
    read_json(ss, result);
    return result;
}

std::vector<double> Calculate_Params(std::vector<double> close)
{
    double price = close[close.size() - 1];

    auto long_term_mean = [&](std::vector<double> x)
    {
        double mu = 0;
        for(auto & mean : x){
            mu += mean;
        }
        mu /= (double) x.size();
        return mu;
    };

    auto long_term_variance = [&](std::vector<double> x)
    {
        double N = (double) x.size();
        double theta = 0;
        double mu = 0;
        for(auto & mean : x){
            mu += mean;
        }
        mu /= N;
        for(auto & value : x){
            theta += pow(value - mu, 2.0);
        }
        theta /= (N - 1);
        return theta;
    };

    auto vol_of_vol = [&](std::vector<double> vol)
    {
        double mean_vol = long_term_mean(vol);
        double volatility = 0;
        for(auto & choice : vol){
            volatility += pow(choice - mean_vol, 2.0);
        }
        volatility /= ((double) vol.size() - 1.0);
        return volatility;
    };

    std::vector<double> result, ror;
    for(int i = 1; i < close.size(); ++i){
        ror.push_back(close[i]/close[i-1] - 1.0);
    }

    double drift = long_term_mean(ror);
    
    int window = 50;
    int NX = ror.size() - window;
    double Theta = 0;
    std::vector<double> mean_variance;

    for(int i = window; i < ror.size(); ++i){
        double tht = long_term_variance({ror.begin() + (i - window), ror.begin() + i});
        mean_variance.push_back(tht);
        Theta += tht;
    }

    Theta /= (double) NX;
    double MU = long_term_mean(mean_variance);

    double top_sum = 0;
    double bot_sum = 0;

    for(int i = 1; i < mean_variance.size(); ++i){
        top_sum += (mean_variance[i] - MU)*(mean_variance[i - 1] - MU);
        bot_sum += pow(mean_variance[i] - MU, 2);
    }

    double p = top_sum/bot_sum;
    double k = -log(p)/1.0;
    
    double sigma = vol_of_vol(mean_variance);
    result.push_back(sqrt(k));
    result.push_back(sqrt(Theta));
    result.push_back(sqrt(sigma));
    result.push_back(price);
    result.push_back(drift);

    return result;
}

double dWT()
{
    int num = 10;
    double dwt = ((rand() % (2*num + 1)) - num)/100.0;
    return dwt;
}

int main()
{
    srand(time(NULL));

    std::string ticker = "MSFT";
    std::string data = GET(ticker);

    std::vector<double> prices;

    Cyclone(JSON(data), std::ref(prices));
    
    std::vector<double> params = Calculate_Params(prices);

    int paths = 100;
    int sims = 1500;

    double t = 22.0/252.0;
    double dt = t / (double) sims;

    double S0, v0, v=0.1;
    
    double final_stock_price = 0;

    for(int i = 0; i < paths; ++i){
        S0 = params[3];
        v0 = v;
        for(int j = 0; j < sims; ++j){
            double dwt = dWT();
            double dwx = -0.5*dwt;
            v0 += params[0]*(params[1] - v0)*dt + params[2]*sqrt(v0)*dwx;
            S0 += params[4]*S0*dt + sqrt(v0)*S0*dwt;
        }
        final_stock_price += S0;
    }

    final_stock_price /= (double) paths;
    
    std::cout << "Your initial stock price is: $" << params[3] << std::endl;
    std::cout << "Your simulated stock price is: $" << final_stock_price << std::endl;

    return 0;
}
