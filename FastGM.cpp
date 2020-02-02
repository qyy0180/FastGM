/*
The code can be further optimized for speed with bit packing and faster random generators.
*/


#include <time.h>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <random>
#include <map>
#include <set>
using namespace std;

vector<default_random_engine> generator_vector;
double threshold = 0.1, v_sum = 0;
int counter = 0;

struct NextBall{
    double x_i;
    int z_i;
    int m_i;
    int c;
};

vector<NextBall> BallVec;
vector<unsigned int> PermutationVec;

void Initial(const vector<double> data, const uint32_t k){
    int n = data.size();
    v_sum = 0;
    NextBall Ball;
    BallVec.resize(n);
    PermutationVec.resize(n*k);

    for(uint32_t i=0;i<n;i++) {
        v_sum = v_sum + data[i];
        Ball.x_i = 0;
        Ball.z_i = 0;
        Ball.m_i = k;
        BallVec[i] = Ball;
        for(int j=0;j<k;j++) PermutationVec[i*k+j] = j;
    }
}

void GetNextBalls(int i, int k){
    default_random_engine generator = generator_vector[i];
    double Q_k = k * 0.2;
    uniform_int_distribution<int> distribution1(1, RAND_MAX-1);
    int r = distribution1(generator);
    double u = (r) / double(RAND_MAX);
    double x=0, z=0;
    int j=0;
    if(BallVec[i].m_i > Q_k){
        x = - log(u);
        j = r % k;
        z = 1;
        counter ++;
    }
    else{
        z =  floor(log(u)/log(1- double(BallVec[i].m_i)/k)) + 1;
        gamma_distribution <double> distribution2(z,1.0);
        x = distribution2(generator);
        j = r % BallVec[i].m_i;
    }
    if(j < BallVec[i].m_i){
        int temp = PermutationVec[i * k + BallVec[i].m_i-1];
        PermutationVec[i * k + BallVec[i].m_i-1] = PermutationVec[i * k + j];
        PermutationVec[i * k + j] = temp;
        BallVec[i].m_i -= 1;
    }

    BallVec[i].x_i += x;
    BallVec[i].z_i += z;
    BallVec[i].c = PermutationVec[i*k + BallVec[i].m_i];
    generator_vector[i] = generator;
}

vector<int> FIPS(const vector<double> data, const uint32_t k, int seed){
    uint32_t delta = 0, k_star = k;
    vector<double> yVec(k, -1);
    vector<int> Sketch(k, -1);

    while (k_star != 0) {
        delta = delta + k;
        for (int i=0;i<data.size();i++){
            double v_i = data[i];
            if(BallVec[i].m_i == k) {
                default_random_engine generator((i << 16) + seed);
                generator_vector.push_back(generator);
            }
            if(v_i == 0) continue;
            uint64_t r_i = delta * v_i / v_sum;

            while(BallVec[i].m_i > 0 and BallVec[i].z_i <= r_i){
                GetNextBalls(i, k);
                double b_i = BallVec[i].x_i /(k*v_i);
                int c = BallVec[i].c;
                if (yVec[c] < 0 ){
                    yVec[c] = b_i;
                    k_star = k_star - 1;
                    Sketch[c] = i;
                }
                else if (b_i < yVec[c]){
                    yVec[c] = b_i;
                    Sketch[c] = i;
                }
            }
        }
    }

    int j_star=0;
    double max = yVec[0];
    for (int i = 0; i < yVec.size(); i++) {
        if (max < yVec[i]){
            max = yVec[i];
            j_star = i;
        }
    }

    //More updates, but faster
    for(int i=0; i!=data.size();i++){
        double v_i = data[i];

        while(BallVec[i].m_i > 0){
            GetNextBalls(i, k);
            int c = BallVec[i].c;

            double b_i = BallVec[i].x_i /(k*v_i);
            if (yVec[j_star] < b_i or BallVec[i].m_i == 0) break;
            else if (b_i < yVec[c]){
                yVec[c] = b_i;
                Sketch[c] = i;
                if(i == j_star){
                    int max = yVec[0];
                    for (int j_i = 0; j_i < yVec.size(); j_i++) {
                        if (max < yVec[j_i]){
                            max = yVec[j_i];
                            j_star = j_i;
                        }
                    }
                }
            }
        }
    }
    generator_vector.clear();
    return Sketch;
}


vector<int> BBM_Mix(const vector<double> data, const uint32_t k, int seed){
    double Q_k = k * threshold;

    vector<vector<double>> yVec;
    vector<double> Vec_temp(k, -1);
    for (int i=0; i<data.size(); i++){
        yVec.push_back(Vec_temp);
    }
    vector<int> Sketch(k, -1);

    for(int i=0; i<data.size(); i++){
        srand((i<<8) + seed);
        default_random_engine generator(((i<<16))+seed);

        vector<int> pi;
        for(int l=0; l<k; l++) pi.push_back(l);

        int m_i = k, z=0;
        double x = 0;
        while(m_i){
            double u = rand() / double(RAND_MAX);
            int j = 0;
            if(m_i > Q_k){
                j = rand() % k;
                x += -log(u) / k / data[i];
            }
            else {
                z = floor(log(u) / log(1-double(m_i)/k)) + 1;
                gamma_distribution <double> distribution (z,1.0);
                x += distribution (generator) / k / data[i];
                j = rand() % m_i;
            }
            if(j < m_i) {
                int temp = pi[j];
                pi[j] = pi[m_i - 1];
                pi[m_i - 1] = temp;
                m_i -= 1;
            }
            yVec[i][pi[m_i]] = x;
        }
    }

    for(int j=0; j<k; j++){
        double min_value = RAND_MAX;
        for(int i=0; i<data.size(); i++){
            if(min_value > yVec[i][j]){
                min_value = yVec[i][j];
                Sketch[j] = i;
            }
        }
    }
    return Sketch;
}

vector<int> BBM_Permutation(const vector<double> data, const uint32_t k, int seed){
    vector<vector<double>> yVec;
    vector<double> Vec_temp(k, -1);
    for (int i=0; i<data.size(); i++){
        yVec.push_back(Vec_temp);
    }
    vector<int> Sketch(k, -1);

    for(int i=0; i<data.size(); i++){
        srand((i<<8) + seed);
        default_random_engine generator(((i<<16))+seed);
        vector<int> pi;
        for(int l=0; l<k; l++) pi.push_back(l);

        int m_i = k, z=0;
        double x = 0;
        while(m_i){
            double u = rand() / double(RAND_MAX);
            int j = 0;
            if(m_i == k) z = 1;
            else z = floor(log(u) / log(1-double(m_i)/k)) + 1;

            gamma_distribution <double> distribution (z,1.0);
            x += distribution (generator) / k / data[i];
            j = rand() % m_i;

            int temp = pi[j];
            pi[j] = pi[m_i-1];
            pi[m_i-1] = temp;
            m_i -= 1;
            yVec[i][pi[m_i]] = x;
        }
    }

    for(int j=0; j<k; j++){
        double min_value = RAND_MAX;
        for(int i=0; i<data.size(); i++){
            if(min_value > yVec[i][j]){
                min_value = yVec[i][j];
                Sketch[j] = i;
            }
        }
    }
    return Sketch;
}

vector<int> BBM_Hash(const vector<double> data, const uint32_t k, int seed){
    vector<vector<double>> yVec;
    vector<double> Vec_temp(k, -1);
    for (int i=0; i<data.size(); i++){
        yVec.push_back(Vec_temp);
    }
    vector<int> Sketch(k, -1);

    for(int i=0; i<data.size(); i++){
        srand((i<<8) + seed);
        int m_i = k;
        double x = 0;
        while(m_i){
            int r = rand();

            double u = r / double(RAND_MAX);
            int j = r % k;
            x += -log(u) / data[i];
            if(yVec[i][j] == -1) {
                yVec[i][j] = x;
                m_i -= 1;
            }
        }
    }

    for(int j=0; j<k; j++){
        double min_value = RAND_MAX;
        for(int i=0; i<data.size(); i++){
            if(min_value > yVec[i][j]){
                min_value = yVec[i][j];
                Sketch[j] = i;
            }
        }
    }
    return Sketch;
}

vector<int> BBM_Basic(const vector<double> data, const uint32_t k, int seed){
    vector<double> yVec(k, -1);
    vector<int> Sketch(k, -1);
    for(int i=0; i<data.size(); i++){
        srand((i<<8) + seed);
        for (int j = 0; j < k; j++) {
            double u = rand() / double(RAND_MAX);
            double exp = -log(u) / data[i];
            if(yVec[j] == -1 or exp < yVec[j]) {
                yVec[j] = exp;
                Sketch[j] = i;
            }
        }
    }
    return Sketch;
}

int main(){
    const uint32_t n = 1000;
    const uint32_t k = 4000;
    const uint32_t round = 1000;
    srand(clock());
    vector<double> TestData;
    for (int i = 0; i<n;i++){
        TestData.emplace_back(rand()/double(RAND_MAX));
    }
    time_t start = clock();

    for (int i = 0; i < round; i++) {
        int seed = i+1;
        Initial(TestData, k);
        vector<int> SketchA = FIPS(TestData, k, seed);

        //vector<int> SketchA = BBM_Mix(TestData, k, seed);

        //vector<int> SketchA = BBM_Permutation(TestData, k, seed);

        //vector<int> SketchA = BBM_Hash(TestData, k, seed);

        //vector<int> SketchA = BBM_Basic(TestData, k, seed);
    }
    time_t end = clock();
    cout << counter / round << endl;
    cout << "Average Time: " << double(end - start) / CLOCKS_PER_SEC / round << endl;
    return 0;
}
