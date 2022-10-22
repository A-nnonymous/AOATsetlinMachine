#include "Habitat.h"
#include "TsetlinMachine.h"
#include "io.h"
#include <climits>

using std::vector;

struct argAndClauses
{
    double                              value;   // Represent precision of this model in testset.
    double                              s;
    int                                 T;
    vector<TsetlinMachine::Clause>      bestNegativeClauses;
    vector<TsetlinMachine::Clause>      bestPositiveClauses;

    argAndClauses()
    {
        value = -(__DBL_MAX__);
    }
    argAndClauses(  double valueIn, double sIn,
                    int TIn,
                    vector<TsetlinMachine::Clause>      bestNeg,
                    vector<TsetlinMachine::Clause>      bestPos)
    {
        value = valueIn;
        s = sIn;
        T = TIn;
        bestNegativeClauses = bestNeg;
        bestPositiveClauses = bestPos;
    }
    
};

struct tsetlinArgs
{
    int numInputs;
    int clausesPerOutput;
    int numOutputs;
    int epochMax;
    vector<double> vars;            // s and T are choosed to be the target.
    std::mt19937 rng;
    

    tsetlinArgs(){}
    tsetlinArgs(int ni, int cpo, int no, double s, int T, int em, std::mt19937 &rng)
    {
        vars.resize(2,0);
        numInputs = ni;
        clausesPerOutput = cpo;
        numOutputs = no;
        vars[0] = s;
        vars[1] = T;
        epochMax = em;
        this->rng = rng;
    }
};


argAndClauses siRNAdemo(tsetlinArgs funcArgs)
{
    std::mt19937                    rng(std::random_device{}());
    int                             train_data_size = 1229;
    int                             test_data_size = 139;
    int                             input_features = funcArgs.numInputs;
    int                             input_bit_per_feature = 4;
    int                             output_size = funcArgs.numOutputs;
    int                             input_size = input_bit_per_feature * input_features;
    int                             clausePerOutput = funcArgs.clausesPerOutput;

    std::vector<std::vector<int>>   train_seqs(train_data_size, std::vector<int>(input_size, 0));
    std::vector<std::vector<int>>   train_scores(train_data_size, std::vector(output_size, 0));
    std::vector<std::vector<int>>   test_seqs(test_data_size, std::vector<int>(input_size, 0));
    std::vector<std::vector<int>>   test_scores(test_data_size, std::vector(output_size, 0));
    parse_huesken_seqs("/home/data/siRNA/e2s/e2s_training_seq.csv", train_seqs);
    parse_huesken_scores("/home/data/siRNA/e2s/e2s_training_efficiency.csv", train_scores);
    parse_huesken_seqs("/home/data/siRNA/e2s/e2s_test_seq.csv", test_seqs);
    parse_huesken_scores("/home/data/siRNA/e2s/e2s_test_efficiency.csv", test_scores);
    std::vector<TsetlinMachine::Clause> bestPositiveClauses;
    std::vector<TsetlinMachine::Clause> bestNegativeClauses;
    std::string outputpath = "/home/output/";
    

    TsetlinMachine  tm( input_size,clausePerOutput,output_size,
                        funcArgs.vars[0],(int)funcArgs.vars[1],
                        rng);
    
    float                           precision; 
    int                             trainset_correct, testset_correct;
    float                           best_test_accuracy = 0;
    std::vector<std::vector<int>>   this_test(test_data_size, std::vector<int>(4,0));
    std::vector<std::vector<int>>   best_test(test_data_size, std::vector<int>(4,0));
    std::vector<int>                hardMaxedResult; hardMaxedResult.resize(4,0);

    for (int epoch = 0; epoch < funcArgs.epochMax; epoch++) {

        trainset_correct = 0; testset_correct = 0;

        for (int trainIdx = 0; trainIdx < train_data_size; trainIdx++)
        {
            TsetlinMachine::Prediction    train_predict = tm.predict(train_seqs[trainIdx]);
            hardMaxedResult = hardMax(train_predict);
            if(hardMaxedResult == train_scores[trainIdx])trainset_correct++;
            tm.learn(train_scores[trainIdx]);
        }
        for(int testIdx = 0; testIdx < test_data_size; testIdx++)
        {
            TsetlinMachine::Prediction    test_predict = tm.predict(test_seqs[testIdx]);
            hardMaxedResult = hardMax(test_predict);
            this_test[testIdx] = hardMaxedResult;
            if(hardMaxedResult == test_scores[testIdx])testset_correct++;
        }

        float this_test_accuracy = testset_correct/(float)test_data_size;
        if(this_test_accuracy > best_test_accuracy)
        {
            best_test=this_test;
            best_test_accuracy = this_test_accuracy;
            tm.getPositiveClauses(&bestPositiveClauses);
            tm.getNegativeClauses(&bestNegativeClauses);
        }
    }
    
    //modelOutput(tm,bestPositiveClauses,bestNegativeClauses,best_test_accuracy, outputpath);
    auto result = argAndClauses(best_test_accuracy,
                                funcArgs.vars[0],funcArgs.vars[1],
                                bestNegativeClauses,
                                bestPositiveClauses);
    return result;
}

int main(int argc, char const *argv[])
{
    // Tsetlin Machine common arguments.
    std::mt19937    rng(std::random_device{}());
    int             inputBitNum = 84;
    int             clausePerOutput = 500;
    int             outputBitNum = 4;
    int             epochNum = 100;
    tsetlinArgs     funcArgs(inputBitNum,clausePerOutput,outputBitNum,0.0,0,epochNum,rng);

    // RSA algorithm arguments;
    int             N = 94;     // Number of individual optimizer.
    int             dimNum = 2;
    int             maxIter = 100;
    double          alpha = 0.1;
    double          beta = 0.005;

    vector<double> mins{0.0, 0.5 * clausePerOutput};
    vector<double> maxes{40.0, 2.0 * clausePerOutput};
    auto limits = Predator<argAndClauses,tsetlinArgs,double>::rangeLimits(mins,maxes);
    auto searchArgs = Predator<argAndClauses,tsetlinArgs,double>::searchArgs(dimNum,maxIter,alpha,beta,limits);
    
    auto RSA = Habitat<argAndClauses,tsetlinArgs,double>(N, siRNAdemo,funcArgs,searchArgs);
    argAndClauses result;
    result = RSA.optimize();
    modelOutput(result.bestPositiveClauses,result.bestNegativeClauses,result.value,"/home/output/");    // Last argument is up to you.

    return 0;
}
