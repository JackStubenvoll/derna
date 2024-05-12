 #include <iostream>
#include <chrono>

#include "Nussinov.h"
#include "NussinovAlgorithm.h"
#include "ZukerAlgorithm.h"
#include "Zuker.h"
#include "default.h"
#include "vector"
#include <string>
#include <tuple>
#include "utils.h"
#include "params/constants.h"

using namespace std;

int main(int argc, char *argv[]) {

    int n;//len of protein
    string input = "../data/uniprotSeq/P15421.fasta";
    string output = "output.txt";
    string rna_file,swipe_output;
    string codon_file = {};
    string param_path = {};
    int model = 1, mode = 1, codonisolate = 0, gen5prime = 0;
    double fphpenergy = -20;
    int fphpposition = 10;
    double target_gc_content = 0.60;
    double incr = inf, lambda = inf, threshold = 0.0025, threshold2 = 0.00075;
    int g = inf;

    if (argc < 2) {
        help();
    }

    try {
        size_t i = 1;
        while ((int)i+1 <= argc) {
            string param = argv[i];
            if (argv[i][0] == '-') {
                switch (argv[i][1]) {
                    case 'i':
                        input = argv[i+1];
                        break;
                    case 'o':
                        output = argv[i+1];
                        break;
                    case 'm':
                        model = std::stoi(argv[i+1]);
                        break;
                    case 's':
                        mode = std::stoi(argv[i+1]);
                        break;
                    case 'g':
                        g = std::stoi(argv[i+1]);
                        break;
                    case 'l':
                        lambda = std::stod(argv[i+1]);
                        break;
                    case 'a':
                        incr = std::stod(argv[i+1]);
                        break;
                    case 'r':
                        rna_file = argv[i+1];
                        break;
                    case 'O':
                        swipe_output = argv[i+1];
                        break;
                    case 'c':
                        codon_file = argv[i+1];
                        break;
                    case 'd':
                        param_path = argv[i+1];
                        break;
                    case 't':
                        threshold = stod(argv[i+1]);
                        break;
                    case 'p':
                        threshold2 = stod(argv[i+1]);
                        break;
                    case 'b':
                        // controls whether or not to isolate first ten codons
                        // only functionality is currently for zuker + mfe_cai
                        codonisolate = std::stoi(argv[i+1]);
                        break;
                    case 'f':
                        //controls whether or not to prepend 5' UTR sequence
                        gen5prime = std::stoi(argv[i+1]);
                        //optional parameters to select 5' UTR sequence
                        if (i + 3 >= argc) {
                            break;
                        }
                        if (!isalpha(argv[i+2][1])) {
                            fphpenergy = stoi(argv[i+2]);
                            fphpposition = stoi(argv[i+3]);
                            if (i + 4 >= argc) {
                                i += 2;
                                break;
                            }
                            if (!isalpha(argv[i+4][1])) {
                                target_gc_content = stod(argv[i+4]);
                                i+=3;
                                break;
                            }
                            i+=2;
                            break;
                        }
                        break;
                    default:
                        help();
                        return(0);
                }
            }
            i += 2;
        }

    } catch (const std::exception& e) {
        std::cout << "Exception!" << std::endl;
        help();
        return -1;
    }
    
    bool nussinov = false;
    bool zuker = false;
    bool test = false;
    switch (model) {
        case 0:
            nussinov = true;
            break;
        case 1:
            zuker = true;
            break;
        case -1:
            test = true;
            break;
        default:
            throw invalid_argument("Invalid Input for Model");
    }

    if (output.empty()) throw invalid_argument("Output File Needed");
    ofstream fout(output);
    scale_params(codon_file, param_path); //"../python/pfizer_codon_usage.csv"


    if (test) {
        if (rna_file.empty()) throw invalid_argument("RNA Input File Needed in Test Mode");
        vector<int> protein = read_fasta(input, fout);

        vector<int> rna = read_rna(rna_file);
        string bp(rna.size(), '.');

        double cai = getCAI(rna, protein);
        double CAI = evaluate_CAI(rna, protein);
        double MFE = evaluate_MFE(rna, bp);

        fout << "secondary structure: " << bp << endl;
        fout << "eval MFE: " << MFE/100 << endl;
        fout << "eval CAI: " << cai << endl;
        fout << "eval standard CAI: " << CAI << endl;
        return 0;
    }

    bool mfe = false;
    bool mfe_cai = false;
    bool lambda_swipe = false;
    bool lambda_swipe2 = false;
    switch (mode) {
        case 1:
            mfe = true;
            break;
        case 2:
            mfe_cai = true;
            break;
        case 3:
            lambda_swipe = true;
            break;
        case 4:
            lambda_swipe2 = true;
            break;
        default:
            throw invalid_argument("Invalid Input for Mode");
    }

    if (gen5prime == 2) {
        cout << "before: " << target_gc_content << endl;
        string temp = generate_Five_Prime(fphpenergy, fphpposition, target_gc_content);
        fout << temp << endl;
        return 0;
    }

    if (input.empty()) throw invalid_argument("Input File Needed");

    vector<int> protein = read_fasta(input, fout);
    n = int(protein.size());
    fout << endl;
    double n_res = 0;
    string rna;

    if (nussinov && mfe) {
        if (g == inf) throw invalid_argument("Invalid Value of g");
        Nussinov N = Nussinov(protein, n, g);
        tuple<double, string, string> temp = N.nussinov(fout);
        n_res = get<0>(temp);
        rna = get<1>(temp);
        int bp = evaluate_BP_N(rna,g);
        fout << "nussinov bp count: " << bp << endl;
    }

    if (nussinov && mfe_cai) {
        if (g == inf) throw invalid_argument("Invalid Value of g");
        Nussinov N = Nussinov(protein, n, g);
        tuple<double, string> temp = N.nussinov_CAI(lambda, fout);
        n_res = get<0>(temp);
        rna = get<1>(temp);
        int bp = evaluate_BP_N(rna,g);
        int type = 0;
        double CAI = evaluate_CAI_N(rna,protein,type);
        fout << "lambda: " << lambda << endl;
        fout << "integrated energy: " << n_res << endl;
        fout << "CAI: " << CAI << endl;
        fout << "nussinov: " << bp << endl;
    }

    if (nussinov && lambda_swipe) {
        Nussinov N = Nussinov(protein, n, g);
        N.lambda_swipe(incr,fout, swipe_output);
    }


    if (zuker && mfe) {
        auto start = chrono::high_resolution_clock::now();
        Zuker Z = Zuker(n,mode,protein);
        Z.calculate_Z(fout);
        Z.traceback_B();
        auto end = chrono::high_resolution_clock::now();
        long time_take = chrono::duration_cast<chrono::seconds>(end - start).count();
        fout << "Time taken : " << time_take;
        fout << "sec" << endl;
        string zuker_bp(3*n,'.'), zuker_rna(3*n,'.'), zuker_rna_X(3*n,'.');// zuker_bp2(3*n,'.'),zuker_bp1(3*n,'.');
        vector<string> rna_array(n, zuker_rna);
        Z.get_bp(zuker_bp);
        Z.get_rna(zuker_rna);
        Z.get_rna_X(zuker_rna_X);
        double cai = evaluate_CAI(zuker_rna, protein, 0);


        fout << "zuker bp:" << zuker_bp << ", size: " << zuker_bp.size() << endl;
        fout << "zuker rna:" << zuker_rna_X << ", size: " << zuker_rna.size() << endl;
        fout << "zuker rna:" << zuker_rna << ", size: " << zuker_rna.size() << endl;
        fout << "zuker cai: " << cai << endl;
        fout << "other rna: " << endl;
    }

    if (zuker && mfe_cai) {
        if (lambda == inf) throw invalid_argument("Invalid Value of lambda");
        if (!codonisolate) {
            Zuker Z = Zuker(n,mode,protein);
            fout << "lambda: " << lambda << endl;
            double energy_cai = Z.calculate_CAI_O(fout, lambda);

            Z.traceback_B2(lambda);
            string zuker_cai_rna(3*n,'.'), zuker_cai_bp(3*n,'.');
            string zuker_cai_rna_X(3*n, '.');
            Z.get_rna_X(zuker_cai_rna_X);
            Z.get_rna_cai(zuker_cai_rna);
            Z.get_bp(zuker_cai_bp);


            int type = 0;

            double CAI_s = evaluate_CAI(zuker_cai_rna,protein,type);
            double CAI = evaluate_CAI(zuker_cai_rna,protein,1);
            double MFE = evaluate_MFE(zuker_cai_rna);

            cout << "lambda: " << lambda << ",O: " << energy_cai << ",cai: " << CAI << ",cai_s: " << CAI_s << ",mfe: " << MFE << ",combined: " << lambda*MFE+(lambda-1)*CAI << endl;
            fout << "zuker cai bp: " << zuker_cai_bp << ",size: " << zuker_cai_bp.size() << endl;
            fout << "zuker rna: " << zuker_cai_rna_X << ".size: " << zuker_cai_rna.size() << endl;
            fout << "zuker cai rna: " << zuker_cai_rna << ".size: " << zuker_cai_rna.size() << endl;

            fout << "Codon Adaptation Index: " << CAI_s << endl;
            fout << "Minimum Free Energy: " << MFE/100 << endl;
        } else if (codonisolate) {
            // for (auto & i : protein) {
            //     cout << i << " "; 
            // }
            vector<int> firstten(10);
            for (int i = 0; i < 10; i++) {
                firstten[i] = protein[i];
            }
            int newn = n-10;
            vector<int> rest(newn);
            for (int i = 10; i < n; i++) {
                rest[i-10] = protein[i];
            }
            // copy(protein.begin(), protein.begin() + 10, firstten.begin());
            // copy(protein.begin() + 10, protein.end(), rest.begin());
            // std::cout << "firstten test" << std::endl;
            // for (int i = 0; i < 10; i++) {
            //     std::cout << firstten[i] << " ";
            // }
            // std::cout << std::endl;
            // for (int i = 0; i < n - 10; i++) {
            //     std::cout << index_aa(rest[i]) << " ";
            // }
            // for (auto & i : rest) {
            //     cout << i << " "; 
            // }
            fout << "first ten amino acids: ";
            for (auto & i : firstten) {
                fout << index_aa(i);
            }
            fout << endl << "rest of amino acid sequence: ";
            for (auto & i : rest) {
                fout << index_aa(i); 
            }
            fout << endl;  
            // Zuker Z = Zuker(n,mode,protein);
            Zuker Z1 = Zuker(10, mode, firstten);
            Zuker Z2 = Zuker(newn, mode, rest);
            
            fout << "lambda: " << lambda << endl;


            double energy_cai_1 = Z1.calculate_CAI_O(fout, lambda);

            Z1.traceback_B2(lambda);
            string zuker_cai_rna_1(3*10,'.'), zuker_cai_bp_1(3*10,'.');
            string zuker_cai_rna_X_1(3*10, '.');
            Z1.get_rna_X(zuker_cai_rna_X_1);
            Z1.get_rna_cai(zuker_cai_rna_1);
            Z1.get_bp(zuker_cai_bp_1);

            int type = 0;


            double energy_cai_2 = Z2.calculate_CAI_O(fout, lambda);

            Z2.traceback_B2(lambda);
            string zuker_cai_rna_2(3*newn,'.'), zuker_cai_bp_2(3*newn,'.');
            string zuker_cai_rna_X_2(3*newn, '.');
            Z2.get_rna_X(zuker_cai_rna_X_2);
            Z2.get_rna_cai(zuker_cai_rna_2);
            Z2.get_bp(zuker_cai_bp_2);

            string zuker_cai_rna_c = zuker_cai_rna_1 + zuker_cai_rna_2;
            string zuker_cai_rna_X_c = zuker_cai_rna_X_1 + zuker_cai_rna_X_2;
            string zuker_cai_bp_c = zuker_cai_bp_1 + zuker_cai_bp_2;

            double CAI_s_1 = evaluate_CAI(zuker_cai_rna_1,firstten,type);
            double CAI_1 = evaluate_CAI(zuker_cai_rna_1,firstten,1);
            double MFE_1 = evaluate_MFE(zuker_cai_rna_1);

            double CAI_s_2 = evaluate_CAI(zuker_cai_rna_2,rest,type);
            double CAI_2 = evaluate_CAI(zuker_cai_rna_2,rest,1);
            double MFE_2 = evaluate_MFE(zuker_cai_rna_2);

            double CAI_s_c = evaluate_CAI(zuker_cai_rna_c,protein,type);
            double CAI_c = evaluate_CAI(zuker_cai_rna_c,protein,1);
            double MFE_c = evaluate_MFE(zuker_cai_rna_c);

            cout << "lambda: " << lambda << endl; 
            cout << ",O_1: " << energy_cai_1 << ",cai_1: " << CAI_1 << ",cai_s_1: " << CAI_s_1 << ",mfe_1: " << MFE_1 << ",combined: " << lambda*MFE_1+(lambda-1)*CAI_1 << endl;
            cout << ",O_2: " << energy_cai_2 << ",cai_2: " << CAI_2 << ",cai_s_2: " << CAI_s_2 << ",mfe_2: " << MFE_2 << ",combined: " << lambda*MFE_2+(lambda-1)*CAI_2 << endl;
            cout << ",O_3: (sum of above): " << energy_cai_1 + energy_cai_2 << ",cai_c: " << CAI_c << ",cai_s_c: " << CAI_s_c << ",mfe_c: " << MFE_c << ",combined: " << lambda*MFE_c+(lambda-1)*CAI_c << endl;
            fout << "zuker cai bp 1: " << zuker_cai_bp_1 << ",size: " << zuker_cai_bp_1.size() << endl;
            fout << "zuker cai bp 2: " << zuker_cai_bp_2 << ",size: " << zuker_cai_bp_2.size() << endl;
            fout << "zuker cai concatenated: " << zuker_cai_bp_c << ",size: " << zuker_cai_bp_c.size() << endl; 
            fout << "zuker rna 1: " << zuker_cai_rna_X_1 << ".size: " << zuker_cai_rna_1.size() << endl;
            fout << "zuker rna 2: " << zuker_cai_rna_X_2 << ".size: " << zuker_cai_rna_2.size() << endl;
            fout << "zuker rna concatenated: " << zuker_cai_rna_X_c << ",size: " << zuker_cai_rna_X_c.size() << endl;
            fout << "zuker cai rna 1: " << zuker_cai_rna_1 << ".size: " << zuker_cai_rna_1.size() << endl;
            fout << "zuker cai rna 1: " << zuker_cai_rna_2 << ".size: " << zuker_cai_rna_2.size() << endl;
            fout << "zuker rna concatenated: " << zuker_cai_rna_c << ",size: " << zuker_cai_rna_c.size() << endl;
            

            fout << "Codon Adaptation Index 1: " << CAI_s_1 << endl;
            fout << "Codon Adaptation Index 2: " << CAI_s_2 << endl;
            fout << "Codon Adaptation Index C: " << CAI_s_c << endl;
            fout << "Minimum Free Energy 1: " << MFE_1/100 << endl;
            fout << "Minimum Free Energy 2: " << MFE_2/100 << endl;
            fout << "Minimum Free Energy C: " << MFE_c/100 << endl;
            rna = zuker_cai_rna_c;
        }
        
    }

    if (zuker && lambda_swipe) {
        Zuker Z = Zuker(n,mode,protein);
        Z.lambda_swipe_2(threshold,threshold2, fout,swipe_output);
    }

    if (zuker && lambda_swipe2) {
        if (lambda == inf) throw invalid_argument("Invalid Value of lambda");
        if (incr == inf) throw invalid_argument("Invalid increment");
        Zuker Z = Zuker(n,mode,protein);
        Z.lambda_swipe(incr,fout,swipe_output);
    }

    if (gen5prime) {
        string fiveprime = generate_Five_Prime(fphpenergy, fphpposition, target_gc_content);
        fout << "five prime seqeunce: " << fiveprime << endl;
        fout << "end sequence: " << fiveprime << rna << endl;

    } 

    return 0;
}





